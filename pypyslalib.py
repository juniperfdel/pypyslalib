import numpy as np


class SLALib:
	# 2 * PI
	D2PI = 6.2831853071795864769252867665590057683943387987502
	DPI = 3.141592653589793238462643e0
	# PI/(12 * 3600) : Seconds to Radians
	DS2R = 7.2722052166430399038487115353692196393452995355905e-5
	# Arcseconds to Radians
	AS2R = 0.484813681109535994e-5

	@staticmethod
	def dmod(A, B):
		return A % B

	@classmethod
	def dranrm(cls, in_val):
		return cls.dmod(in_val, cls.D2PI)

	@classmethod
	def gmst(cls, in_ut1):
		# Julian centuries from fundamental epoch J2000 to this UT
		tu = (in_ut1 - 51544.5) / 36525.0
		return cls.dranrm(cls.dmod(in_ut1, 1.0) * cls.D2PI + (24110.54841 + (8640184.812866 + (0.093104 - 6.2e-6 * tu) * tu) * tu) * cls.DS2R)

	@classmethod
	def hour_angle(cls, in_mjd, in_ra, in_long):
		# Not part of the original library
		return np.rad2deg(cls.dranrm(cls.dranrm(cls.gmst(in_mjd) + in_long) - in_ra))

	@classmethod
	def clyd(cls, input_year, input_month, input_day):
		# +
		# - - - - -
		# C L Y D
		# - - - - -
		# Gregorian calendar to year and day in year (in a Julian calendar
		# aligned to the 20th/21st century Gregorian calendar).
		# Given:
		# IY,IM,ID   i    year, month, day in Gregorian calendar
		# Returned:
		# NY         i    year (re-aligned Julian calendar)
		# ND         i    day in year (1 = January 1st)
		# JSTAT      i    status:
		# 0 = OK
		# 1 = bad year (before -4711)
		# 2 = bad month
		# 3 = bad day (but conversion performed)
		# Notes:
		# 1  This routine exists to support the low-precision routines
		# sla_EARTH, sla_MOON and sla_ECOR.
		# 2  Between 1900 March 1 and 2100 February 28 it returns answers
		# which are consistent with the ordinary Gregorian calendar.
		# Outside this range there will be a discrepancy which increases
		# by one day for every non-leap century year.
		# 3  The essence of the algorithm is first to express the Gregorian
		# date as a Julian Day Number and then to convert this back to
		# a Julian calendar date, with day-in-year instead of month and
		# day.  See 12.92-1 and 12.95-1 in the reference.
		# Reference:  Explanatory Supplement to the Astronomical Almanac,
		# ed P.K.Seidelmann, University Science Books (1992),
		# p604-606.

		return_code = 0
		return_year = 0
		return_day = 0

		# Validate year
		if input_year >= -4711:
			#  Validate month
			if 1 <= input_month <= 12:
				month_lengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

				# Allow for (Gregorian) leap year
				if np.mod(input_year, 4) == 0 and (np.mod(input_year, 100) != 0 or np.mod(input_year, 400) == 0):
					month_lengths[1] = 29

				# Validate day
				if input_day < 1 or input_day > month_lengths[int(input_month) - 1]:
					return_code = 3
				# Perform the conversion
				temp_i = (14 - input_month) / 12
				temp_k = input_year - temp_i
				temp_j = (1461 * (temp_k + 4800)) / 4 + (367 * (input_month - 2 + 12 * temp_i)) / 12 - (3 * ((temp_k + 4900) / 100)) / 4 + input_day - 30660
				temp_k = (temp_j - 1) / 1461
				temp_l = temp_j - 1461 * temp_k
				temp_n = (temp_l - 1) / 365 - temp_l / 1461
				temp_j = ((80 * (temp_l - 365 * temp_n + 30)) / 2447) / 11
				temp_i = temp_n + temp_j
				return_day = 59 + temp_l - 365 * temp_i + ((4 - temp_n) / 4) * (1 - temp_j)
				return_year = 4 * temp_k + temp_i - 4716

			# Bad month
			else:
				return_code = 2
		else:
			# Bad year
			return_code = 1

		return return_year, return_day, return_code

	@classmethod
	def dcc2s(cls, in_coord):
		# +
		#      - - - - - -
		#       D C C 2 S
		#      - - - - - -
		#   Cartesian to spherical coordinates
		#   Given:
		#      V     d(3)   x,y,z vector
		#   Returned:
		#      A,B   d      spherical coordinates in radians
		#   The spherical coordinates are longitude (+ve anticlockwise looking
		#   from the +ve latitude pole) and latitude.  The Cartesian coordinates
		#   are right handed, with the x axis at zero longitude and latitude, and
		#   the z axis at the +ve latitude pole.
		#   If V is null, zero A and B are returned.  At either pole, zero A is
		#   returned.
		#   Last revision:   22 July 2004

		x, y, z = in_coord
		r = np.sqrt(x * x + y * y)

		longitude = np.where(r == 0, 0, np.arctan2(y, x))
		latitude = np.where(z == 0, 0, np.arctan2(z, r))

		return longitude, latitude

	@classmethod
	def dcs2c(cls, ra, dec):
		# +
		#      - - - - - -
		#       D C S 2 C
		#      - - - - - -
		#
		#   Spherical coordinates to direction cosines (double precision)
		#
		#   Given:
		#      A,B       d      spherical coordinates in radians
		#                          (RA,Dec), (long,lat) etc.
		#
		#   Returned:
		#      V         d(3)   x,y,z unit vector
		#
		#   The spherical coordinates are longitude (+ve anticlockwise looking
		#   from the +ve latitude pole) and latitude.  The Cartesian coordinates
		#   are right handed, with the x axis at zero longitude and latitude, and
		#   the z axis at the +ve latitude pole.
		#
		#   Last revision:   26 December 2004

		return np.array([np.cos(ra) * np.cos(dec), np.sin(ra) * np.cos(dec), np.sin(dec)])

	@classmethod
	def deuler(cls, order, phi, theta, psi):
		#
		#      - - - - - - -
		#       D E U L E R
		#      - - - - - - -

		#   Form a rotation matrix from the Euler angles - three successive
		#   rotations about specified Cartesian axes (double precision)

		#   Given:
		#     ORDER   c*(*)   specifies about which axes the rotations occur
		#     PHI     d       1st rotation (radians)
		#     THETA   d       2nd rotation (   "   )
		#     PSI     d       3rd rotation (   "   )

		#   Returned:
		#     RMAT    d(3,3)  rotation matrix

		#   A rotation is positive when the reference frame rotates
		#   anticlockwise as seen looking towards the origin from the
		#   positive region of the specified axis.

		#   The characters of ORDER define which axes the three successive
		#   rotations are about.  A typical value is 'ZXZ', indicating that
		#   RMAT is to become the direction cosine matrix corresponding to
		#   rotations of the reference frame through PHI radians about the
		#   old Z-axis, followed by THETA radians about the resulting X-axis,
		#   then PSI radians about the resulting Z-axis.

		#   The axis names can be any of the following, in any order or
		#   combination:  X, Y, Z, uppercase or lowercase, 1, 2, 3.  Normal
		#   axis labelling/numbering conventions apply;  the xyz (=123)
		#   triad is right-handed.  Thus, the 'ZXZ' example given above
		#   could be written 'zxz' or '313' (or even 'ZxZ' or '3xZ').  ORDER
		#   is terminated by length or by the first unrecognized character.

		#   Fewer than three rotations are acceptable, in which case the later
		#   angle arguments are ignored.  If all rotations are zero, the
		#   identity matrix is produced.

		#   Initialize result matrix
		result_mat = np.identity(3)
		#   Look at each character of axis string until finished
		for i_char, axis in enumerate(order):
			# Initialize rotation matrix for the current rotation
			rot_mat = np.identity(3)
			# Pick up the appropriate Euler angle and take sine & cosine
			if i_char == 1:
				angle = phi
			elif i_char == 2:
				angle = theta
			else:
				angle = psi

			ang_sin = np.sin(angle)
			ang_cos = np.cos(angle)

			# Identify the axis
			if axis in ['X', 'x', '1']:
				# Matrix for x-rotation
				rot_mat[1, 1] = ang_cos
				rot_mat[1, 2] = ang_sin
				rot_mat[2, 1] = -ang_sin
				rot_mat[2, 2] = ang_cos

			elif axis in ['Y', 'y', '2']:
				# Matrix for y-rotation
				rot_mat[0, 0] = ang_cos
				rot_mat[0, 2] = -ang_sin
				rot_mat[2, 0] = ang_sin
				rot_mat[2, 2] = ang_cos

			elif axis in ['Z', 'z', '3']:
				# Matrix for z-rotation
				rot_mat[0, 0] = ang_cos
				rot_mat[0, 1] = ang_sin
				rot_mat[1, 0] = -ang_sin
				rot_mat[1, 1] = ang_cos
			else:
				raise ValueError("Invalid Character received for deuler!")

			result_mat = result_mat @ rot_mat
		return result_mat

	@classmethod
	def djcl(cls, in_mjd):
		# +
		# - - - - -
		# D J C L
		# - - - - -
		# Modified Julian Date to Gregorian year, month, day,
		# and fraction of a day.
		# Given:
		# DJM      dp     modified Julian Date (JD-2400000.5)
		# Returned:
		# IY       int    year
		# IM       int    month
		# ID       int    day
		# FD       dp     fraction of day
		# J        int    status:
		# 0 = OK
		# -1 = unacceptable date (before 4701BC March 1)
		# The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).
		# Last revision:   22 July 2004

		# Check if date is acceptable.

		r_year = r_month = r_day = r_frac_day = 0

		if in_mjd <= -2395520e0 or in_mjd >= 1e9:
			return_code = -1
		else:
			# Separate day and fraction.
			frac_day = np.mod(in_mjd, 1e0)
			if frac_day < 0e0:  frac_day += 1e0
			day_int = np.rint(in_mjd - frac_day)

			# Express day in Gregorian calendar.
			jd = np.rint(day_int) + 2400001

			n4 = 4 * (jd + ((6 * ((4 * jd - 17918) / 146097)) / 4 + 1) / 2 - 37)
			nd10 = 10 * (np.mod(n4 - 237, 1461) / 4) + 5

			r_year = n4 / 1461 - 4712
			r_month = np.mod(nd10 / 306 + 2, 12) + 1
			r_day = np.mod(nd10, 306) / 10 + 1
			r_frac_day = frac_day

			return_code = 0
		return r_year, r_month, r_day, r_frac_day, return_code

	@classmethod
	def dmxv(cls, in_mat, in_vec):
		# +
		# - - - - -
		# D M X V
		# - - - - -
		#
		# Performs the 3-D forward unitary transformation:
		#
		# vector VB = matrix DM * vector VA
		#
		# (double precision)
		#
		# Given:
		# DM       dp(3,3)    matrix
		# VA       dp(3)      vector
		#
		# Returned:
		# VB       dp(3)      result vector
		#
		# To comply with the ANSI Fortran 77 standard, VA and VB must be
		# different arrays.  However, the routine is coded so as to work
		# properly on many platforms even if this rule is violated.
		#
		# Last revision:   26 December 2004
		# -

		return np.dot(in_mat, in_vec)

	@classmethod
	def prebn(cls, bess_epoch_0, bess_epoch_1):
		# +
		# - - - - - -
		# P R E B N
		# - - - - - -
		#
		# Generate the matrix of precession between two epochs,
		# using the old, pre-IAU1976, Bessel-Newcomb model, using
		# Kinoshita's formulation (double precision)
		#
		# Given:
		# BEP0    dp         beginning Besselian epoch
		# BEP1    dp         ending Besselian epoch
		#
		# Returned:
		# RMATP  dp(3,3)    precession matrix
		#
		# The matrix is in the sense   V(BEP1)  =  RMATP * V(BEP0)
		#
		# Reference:
		# Kinoshita, H. (1975) 'Formulas for precession', SAO Special
		# Report No. 364, Smithsonian Institution Astrophysical
		# Observatory, Cambridge, Massachusetts.
		#
		# Called:  sla_DEULER
		# -

		# Interval between basic epoch B1850.0 and beginning epoch in TC
		epoch_diff = (bess_epoch_0 - 1850e0) / 100e0

		# Interval over which precession required, in tropical centuries
		interval_precision = (bess_epoch_1 - bess_epoch_0) / 100e0

		# Euler angles
		tas2_r = interval_precision * cls.AS2R
		w = 2303.5548e0 + (1.39720e0 + 0.000059e0 * epoch_diff) * epoch_diff

		zeta = (w + (0.30242e0 - 0.000269e0 * epoch_diff + 0.017996e0 * interval_precision) * interval_precision) * tas2_r
		z = (w + (1.09478e0 + 0.000387e0 * epoch_diff + 0.018324e0 * interval_precision) * interval_precision) * tas2_r
		theta = (2005.1125e0 + (-0.85294e0 - 0.000365e0 * epoch_diff) * epoch_diff + (-0.42647e0 - 0.000365e0 * epoch_diff - 0.041802e0 * interval_precision) * interval_precision) * tas2_r

		# Rotation matrix
		return cls.deuler('ZYZ', -zeta, theta, -z)

	@classmethod
	def prec(cls, start_epoch, end_epoch):
		# +
		# - - - - -
		# P R E C
		# - - - - -
		#
		# Form the matrix of precession between two epochs (IAU 1976, FK5)
		# (double precision)
		#
		# Given:
		# EP0    dp         beginning epoch
		# EP1    dp         ending epoch
		#
		# Returned:
		# RMATP  dp(3,3)    precession matrix
		#
		# Notes:
		#
		# 1)  The epochs are TDB (loosely ET) Julian epochs.
		#
		# 2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
		#
		# 3)  Though the matrix method itself is rigorous, the precession
		# angles are expressed through canonical polynomials which are
		# valid only for a limited time span.  There are also known
		# errors in the IAU precession rate.  The absolute accuracy
		# of the present formulation is better than 0.1 arcsec from
		# 1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
		# and remains below 3 arcsec for the whole of the period
		# 500BC to 3000AD.  The errors exceed 10 arcsec outside the
		# range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
		# 5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
		# The SLALIB routine sla_PRECL implements a more elaborate
		# model which is suitable for problems spanning several
		# thousand years.
		#
		# References:
		# Lieske,J.H., 1979. Astron.Astrophys.,73,282.
		# equations (6) & (7), p283.
		# Kaplan,G.H., 1981. USNO circular no. 163, pA2.
		#
		# Called:  sla_DEULER
		# -

		# Interval between basic epoch J2000.0 and beginning epoch (JC)
		t0 = (start_epoch - 2000e0) / 100e0

		# Interval over which precession required (JC)
		interval_precision = (end_epoch - start_epoch) / 100e0

		# Euler angles
		tas2_r = interval_precision * cls.AS2R
		w = 2306.2181e0 + (1.39656e0 - 0.000139e0 * t0) * t0

		zeta = (w + ((0.30188e0 - 0.000344e0 * t0) + 0.017998e0 * interval_precision) * interval_precision) * tas2_r
		z = (w + ((1.09468e0 + 0.000066e0 * t0) + 0.018203e0 * interval_precision) * interval_precision) * tas2_r
		theta = ((2004.3109e0 + (-0.85330e0 - 0.000217e0 * t0) * t0) + ((-0.42665e0 - 0.000217e0 * t0) - 0.041833e0 * interval_precision) * interval_precision) * tas2_r

		# Rotation matrix
		return cls.deuler('ZYZ', -zeta, theta, -z)

	@classmethod
	def preces(cls, system_name, start_epoch, end_epoch, input_ra, input_dec):
		# +
		# - - - - - - -
		# P R E C E S
		# - - - - - - -
		#
		# Precession - either FK4 (Bessel-Newcomb, pre IAU 1976) or
		# FK5 (Fricke, post IAU 1976) as required.
		#
		# Given:
		# SYSTEM     char   precession to be applied: 'FK4' or 'FK5'
		# EP0,EP1    dp     starting and ending epoch
		# RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP0
		#
		# Returned:
		# RA,DC      dp     RA,Dec, mean equator & equinox of epoch EP1
		#
		# Called:    sla_DRANRM, sla_PREBN, sla_PREC, sla_DCS2C,
		# sla_DMXV, sla_DCC2S
		#
		# Notes:
		#
		# 1)  Lowercase characters in SYSTEM are acceptable.
		#
		# 2)  The epochs are Besselian if SYSTEM='FK4' and Julian if 'FK5'.
		# For example, to precess coordinates in the old system from
		# equinox 1900.0 to 1950.0 the call would be:
		# CALL sla_PRECES ('FK4', 1900D0, 1950D0, RA, DC)
		#
		# 3)  This routine will NOT correctly convert between the old and
		# the new systems - for example conversion from B1950 to J2000.
		# For these purposes see sla_FK425, sla_FK524, sla_FK45Z and
		# sla_FK54Z.
		#
		# 4)  If an invalid SYSTEM is supplied, values of -99,-99 will
		# be returned for both RA and DC.
		# -

		# Convert to uppercase and validate SYSTEM
		if system_name.upper() not in ['FK4', 'FK5']:
			ra = -99e0
			dec = -99e0
		else:
			# Generate appropriate precession matrix
			if system_name.upper() == 'FK4':
				prec_mat = cls.prebn(start_epoch, end_epoch)
			else:
				prec_mat = cls.prec(start_epoch, end_epoch)

			# Convert RA,Dec to x,y,z
			v1 = cls.dcs2c(input_ra, input_dec)

			# Precess
			v2 = cls.dmxv(prec_mat, v1)

			# Back to RA,Dec
			ra, dec = cls.dcc2s(v2)
			ra = cls.dranrm(ra)

		return ra, dec

	@classmethod
	def precl(cls, epoch_start, epoch_end):
		# +
		# - - - - - -
		# P R E C L
		# - - - - - -
		#
		# Form the matrix of precession between two epochs, using the
		# model of Simon et al (1994), which is suitable for long
		# periods of time.
		#
		# (double precision)
		#
		# Given:
		# EP0    dp         beginning epoch
		# EP1    dp         ending epoch
		#
		# Returned:
		# RMATP  dp(3,3)    precession matrix
		#
		# Notes:
		#
		# 1)  The epochs are TDB Julian epochs.
		#
		# 2)  The matrix is in the sense   V(EP1)  =  RMATP * V(EP0)
		#
		# 3)  The absolute accuracy of the model is limited by the
		# uncertainty in the general precession, about 0.3 arcsec per
		# 1000 years.  The remainder of the formulation provides a
		# precision of 1 mas over the interval from 1000AD to 3000AD,
		# 0.1 arcsec from 1000BC to 5000AD and 1 arcsec from
		# 4000BC to 8000AD.
		#
		# Reference:
		# Simon, J.L. et al., 1994. Astron.Astrophys., 282, 663-683.
		#
		# Called:  sla_DEULER
		# -

		# Interval between basic epoch J2000.0 and beginning epoch (1000JY)
		t0 = (epoch_start - 2000e0) / 1000e0

		# Interval over which precession required (1000JY)
		t = (epoch_end - epoch_start) / 1000e0

		# Euler angles
		tas2_r = t * cls.AS2R
		w = 23060.9097e0 + (139.7459e0 + (-0.0038e0 + (-0.5918e0 + (-0.0037e0 + 0.0007e0 * t0) * t0) * t0) * t0) * t0

		zeta = (w + (30.2226e0 + (-0.2523e0 + (-0.3840e0 + (-0.0014e0 + 0.0007e0 * t0) * t0) * t0) * t0 + (
				18.0183e0 + (-0.1326e0 + (0.0006e0 + 0.0005e0 * t0) * t0) * t0 + (-0.0583e0 + (-0.0001e0 + 0.0007e0 * t0) * t0 + (-0.0285e0 + (-0.0002e0) * t) * t) * t) * t) * t) * tas2_r

		z = (w + (109.5270e0 + (0.2446e0 + (-1.3913e0 + (-0.0134e0 + 0.0026e0 * t0) * t0) * t0) * t0 + (18.2667e0 + (-1.1400e0 + (-0.0173e0 + 0.0044e0 * t0) * t0) * t0 + (
				-0.2821e0 + (-0.0093e0 + 0.0032e0 * t0) * t0 + (-0.0301e0 + 0.0006e0 * t0 - 0.0001e0 * t) * t) * t) * t) * t) * tas2_r

		theta = (20042.0207e0 + (-85.3131e0 + (-0.2111e0 + (0.3642e0 + (0.0008e0 + (-0.0005e0) * t0) * t0) * t0) * t0) * t0 + (
				-42.6566e0 + (-0.2111e0 + (0.5463e0 + (0.0017e0 + (-0.0012e0) * t0) * t0) * t0) * t0 + (-41.8238e0 + (0.0359e0 + (0.0027e0 + (-0.0001e0) * t0) * t0) * t0 + (
				-0.0731e0 + (0.0019e0 + 0.0009e0 * t0) * t0 + (-0.0127e0 + 0.0011e0 * t0 + 0.0004e0 * t) * t) * t) * t) * t) * tas2_r

		# Rotation matrix
		return cls.deuler('ZYZ', -zeta, theta, -z)

	@classmethod
	def dh2e(cls, az, el, phi):
		# +
		# - - - - -
		# D E 2 H
		# - - - - -
		#
		# Horizon to equatorial coordinates:  Az,El to HA,Dec
		#
		# (double precision)
		#
		# Given:
		# AZ      d     azimuth
		# EL      d     elevation
		# PHI     d     observatory latitude
		#
		# Returned:
		# HA      d     hour angle
		# DEC     d     declination
		#
		# Notes:
		#
		# 1)  All the arguments are angles in radians.
		#
		# 2)  The sign convention for azimuth is north zero, east +pi/2.
		#
		# 3)  HA is returned in the range +/-pi.  Declination is returned
		# in the range +/-pi/2.
		#
		# 4)  The latitude is (in principle) geodetic.  In critical
		# applications, corrections for polar motion should be applied.
		#
		# 5)  In some applications it will be important to specify the
		# correct type of elevation in order to produce the required
		# type of HA,Dec.  In particular, it may be important to
		# distinguish between the elevation as affected by refraction,
		# which will yield the "observed" HA,Dec, and the elevation
		# in vacuo, which will yield the "topocentric" HA,Dec.  If the
		# effects of diurnal aberration can be neglected, the
		# topocentric HA,Dec may be used as an approximation to the
		# "apparent" HA,Dec.
		#
		# 6)  No range checking of arguments is done.
		#
		# 7)  In applications which involve many such calculations, rather
		# than calling the present routine it will be more efficient to
		# use inline code, having previously computed fixed terms such
		# as sine and cosine of latitude.
		# -

		# Useful trig functions
		sa = np.sin(az)
		ca = np.cos(az)
		se = np.sin(el)
		ce = np.cos(el)
		sp = np.sin(phi)
		cp = np.cos(phi)

		# HA,Dec as x,y,z
		x = -ca * ce * sp + se * cp
		y = -sa * ce
		z = ca * ce * cp + se * sp

		# To HA,Dec
		r = np.sqrt(x * x + y * y)

		ha = np.where(r == 0, 0, np.arctan2(y, x))
		dec = np.arctan2(z, r)
		return ha, dec

	@classmethod
	def gmsta(cls, date, ut=None):
		# +
		# - - - - - -
		# G M S T A
		# - - - - - -
		#
		# Conversion from Universal Time to Greenwich mean sidereal time,
		# with rounding errors minimized.
		#
		# double precision
		#
		# Given:
		# DATE    d      UT1 date (MJD: integer part of JD-2400000.5))
		# UT      d      UT1 time (fraction of a day)
		#
		# The result is the Greenwich mean sidereal time (double precision,
		# radians, in the range 0 to 2pi).
		#
		# There is no restriction on how the UT is apportioned between the
		# DATE and UT arguments.  Either of the two arguments could, for
		# example, be zero and the entire date+time supplied in the other.
		# However, the routine is designed to deliver maximum accuracy when
		# the DATE argument is a whole number and the UT lies in the range
		# 0 to 1 (or vice versa).
		#
		# The algorithm is based on the IAU 1982 expression (see page S15 of
		# the 1984 Astronomical Almanac).  This is always described as giving
		# the GMST at 0 hours UT1.  In fact, it gives the difference between
		# the GMST and the UT, the steady 4-minutes-per-day drawing-ahead of
		# ST with respect to UT.  When whole days are ignored, the expression
		# happens to equal the GMST at 0 hours UT1 each day.  Note that the
		# factor 1.0027379... does not appear explicitly but in the form of
		# the coefficient 8640184.812866, which is 86400x36525x0.0027379...
		#
		# In this routine, the entire UT1 (the sum of the two arguments DATE
		# and UT) is used directly as the argument for the standard formula.
		# The UT1 is then added, but omitting whole days to conserve accuracy.
		#
		# See also the routine sla_GMST, which accepts the UT as a single
		# argument.  Compared with sla_GMST, the extra numerical precision
		# delivered by the present routine is unlikely to be important in
		# an absolute sense, but may be useful when critically comparing
		# algorithms and in applications where two sidereal times close
		# together are differenced.
		#
		# Called:  sla_DRANRM
		# -

		if ut is None:
			ut = np.mod(date, 1.0)
			date = np.trunc(date)

		# Seconds of time to radians

		# Julian centuries since J2000.
		if date < ut:
			d1 = date
			d2 = ut
		else:
			d1 = ut
			d2 = date

		T = (d1 + (d2 - 51544.5e0)) / 36525e0

		# GMST at this UT1.
		return cls.dranrm(cls.DS2R * (24110.54841e0 + (8640184.812866e0 + (0.093104e0 - 6.2e-6 * T) * T) * T + 86400e0 * (np.mod(d1, 1e0) + np.mod(d2, 1e0))))

	@classmethod
	def eqgal(cls, in_ra, in_dec):
		# +
		# - - - - - -
		# E Q G A L
		# - - - - - -
		#
		# Transformation from J2000.0 equatorial coordinates to
		# IAU 1958 galactic coordinates (double precision)
		#
		# Given:
		# DR,DD       dp       J2000.0 RA,Dec
		#
		# Returned:
		# DL,DB       dp       galactic longitude and latitude L2,B2
		#
		# (all arguments are radians)
		#
		# Called:
		# sla_DCS2C, sla_DMXV, sla_DCC2S, sla_DRANRM, sla_DRANGE
		#
		# Note:
		# The equatorial coordinates are J2000.0.  Use the routine
		# sla_EG50 if conversion from B1950.0 'FK4' coordinates is
		# required.

		# L2,B2 system of galactic coordinates
		#
		# P = 192.25       RA of galactic north pole (mean B1950.0)
		# Q =  62.6        inclination of galactic to mean B1950.0 equator
		# R =  33          longitude of ascending node
		#
		# P,Q,R are degrees
		#
		# Equatorial to galactic rotation matrix (J2000.0), obtained by
		# applying the standard FK4 to FK5 transformation, for zero proper
		# motion in FK5, to the columns of the B1950 equatorial to
		# galactic rotation matrix:
		#

		RMAT = np.array(
			[
				[-0.054875539726e0, -0.873437108010e0, -0.483834985808e0],
				[0.494109453312e0, -0.444829589425e0, 0.746982251810e0],
				[-0.867666135858e0, -0.198076386122e0, 0.455983795705e0]
			]
		)

		# Spherical to Cartesian
		V1 = cls.dcs2c(in_ra, in_dec)

		# Equatorial to galactic
		V2 = cls.dmxv(RMAT, V1)

		# Cartesian to spherical
		r_gal_l, r_gal_b = cls.dcc2s(V2)

		# Express in conventional ranges
		r_gal_l = cls.dranrm(r_gal_l)
		r_gal_b = cls.drange(r_gal_b)

		return r_gal_l, r_gal_b

	@classmethod
	def drange(cls, in_angle):
		rv = np.mod(in_angle, cls.D2PI)
		rv = np.where(np.abs(rv) >= cls.DPI, rv - np.sign(in_angle) * cls.D2PI, rv)
		return rv
