import numpy as np


class SLALib:
    # 2 * PI
    D2PI = 6.2831853071795864769252867665590057683943387987502
    # PI
    DPI = 3.141592653589793238462643
    # PI/(12 * 3600) : Seconds to Radians
    DS2R = 7.2722052166430399038487115353692196393452995355905e-5
    # Arcseconds to Radians
    AS2R = 4.848136811095359935899141e-6
    # Epsilon for float comparision
    TINY = 1e-30
    # Light time for 1 AU (sec)
    CR = 499.004782e0
    # Gravitational radius of the Sun x 2 (2*mu/c**2, AU)
    GR2 = 1.974126e-8
    # B1950
    B1950 = 1949.9997904423e0
    # Degrees to radians
    DD2R = 1.745329251994329576923691e-2
    # Arc seconds in a full circle
    TURNAS = 1296000e0
    # Reference epoch (J2000), MJD
    DJM0 = 51544.5e0
    # Days per Julian century
    DJC = 36525e0
    # Mean sidereal rate (at J2000) in radians per (UT1) second
    SR = 7.292115855306589e-5
    # Earth equatorial radius (metres)
    A0 = 6378140e0
    # Reference spheroid flattening factor and useful function
    SPHF = 1e0 / 298.257e0
    SPHB = (1e0 - SPHF) ** 2
    # Astronomical unit in metres
    AU=1.49597870e11

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
            if frac_day < 0e0:
                frac_day += 1e0
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
    def djcal(cls, ndp, djm):
        """
        - - - - - -
        D J C A L
        - - - - - -

        Modified Julian Date to Gregorian Calendar, expressed
        in a form convenient for formatting messages (namely
        rounded to a specified precision, and with the fields
        stored in a single array)

        Given:
        NDP      i      number of decimal places of days in fraction
        DJM      d      modified Julian Date (JD-2400000.5)

        Returned:
        IYMDF    i(4)   year, month, day, fraction in Gregorian
        calendar
        J        i      status:  nonzero = out of range

        Any date after 4701BC March 1 is accepted.

        NDP should be 4 or less if internal overflows are to be avoided
        on machines which use 32-bit integers.

        The algorithm is adapted from Hatcher 1984 (QJRAS 25, 53-55).

        Last revision:   22 July 2004
        """
        # Validate.
        iymdf = np.array([])
        if (djm <= -2395520e0) or (djm <= -1e9):
            j = -1
        else:
            j = 0
        
            # Denominator of fraction.
            NFD = 10 ** np.maximum(ndp, 0)
            FD = NFD
        
            # Round date and express in units of fraction.
            DF = np.rint(djm * FD)
        
            # Separate day and fraction.
            F = np.mod(DF, FD)
            if F < 0e0:
                F = F + FD
            D = (DF - F) / FD
        
            # Express day in Gregorian calendar.
            JD = np.rint(D) + 2400001
        
            N4 = 4 * (JD + ((2 * ((4 * JD - 17918) / 146097) * 3) / 4 + 1) / 2 - 37)
            ND10 = 10 * (np.mod(N4 - 237, 1461) / 4) + 5
        
            iymdf = np.array([N4 / 1461 - 4712, np.mod(ND10 / 306 + 2, 12) + 1, np.mod(ND10, 306) / 10 + 1, np.rint(F)])
        return j, iymdf

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
    def dimxv(cls, in_mat, in_vec):
        """
        - - - - - -
        D I M X V
        - - - - - -

        Performs the 3-D backward unitary transformation:

        vector VB = (inverse of matrix DM) * vector VA

        (double precision)

        (n.b.  the matrix must be unitary, as this routine assumes that
        the inverse and transpose are identical)

        Given:
        DM       dp(3,3)    matrix
        VA       dp(3)      vector

        Returned:
        VB       dp(3)      result vector
        Inverse of matrix DM * vector VA -> vector VW
        """
        return np.dot(in_mat.T, in_vec)

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

        theta = (20042.0207e0 + (-85.3131e0 + (-0.2111e0 + (0.3642e0 + (0.0008e0 + (-0.0005e0) * t0) * t0) * t0) * t0) * t0 + (-42.6566e0 + (-0.2111e0 + (0.5463e0 + (0.0017e0 + (-0.0012e0) * t0) * t0) * t0) *
                 t0 + (-41.8238e0 + (0.0359e0 + (0.0027e0 + (-0.0001e0) * t0) * t0) * t0 + (-0.0731e0 + (0.0019e0 + 0.0009e0 * t0) * t0 + (-0.0127e0 + 0.0011e0 * t0 + 0.0004e0 * t) * t) * t) * t) * t) * tas2_r

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
        date_lt_ut = date < ut
        d1 = np.where(date_lt_ut, date, ut)
        d2 = np.where(date_lt_ut, ut, date)

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

    @classmethod
    def addet(cls, rm, dm, eq):
        # +
        # - - - - - -
        # A D D E T
        # - - - - - -
        #
        # Add the E-terms (elliptic component of annual aberration)
        # to a pre IAU 1976 mean place to conform to the old
        # catalogue convention (double precision)
        #
        # Given:
        # RM,DM     dp     RA,Dec (radians) without E-terms
        # EQ        dp     Besselian epoch of mean equator and equinox
        #
        # Returned:
        # RC,DC     dp     RA,Dec (radians) with E-terms included
        #
        # Note:
        #
        # Most star positions from pre-1984 optical catalogues (or
        # derived from astrometry using such stars) embody the
        # E-terms.  If it is necessary to convert a formal mean
        # place (for example a pulsar timing position) to one
        # consistent with such a star catalogue, then the RA,Dec
        # should be adjusted using this routine.
        #
        # Reference:
        # Explanatory Supplement to the Astronomical Ephemeris,
        # section 2D, page 48.
        #
        # Depends:  sla_ETRMS, sla_DCS2C, sla_DCC2S, sla_DRANRM, sla_DRANGE

        # E-terms vector
        A = cls.etrms(eq)

        # Spherical to Cartesian
        V = cls.dcs2c(rm, dm)

        # Include the E-terms
        V += A

        # Cartesian to spherical
        rc, dc = cls.dcc2s(V)

        # Bring RA into conventional range
        rc = cls.dranrm(rc)

        return rc, dc

    @classmethod
    def etrms(cls, ep):
        # +
        # - - - - - -
        # E T R M S
        # - - - - - -
        #
        # Compute the E-terms (elliptic component of annual aberration)
        # vector (double precision)
        #
        # Given:
        # EP      dp      Besselian epoch
        #
        # Returned:
        # EV      dp(3)   E-terms as (dx,dy,dz)
        #
        # Note the use of the J2000 aberration constant (20.49552 arcsec).
        # This is a reflection of the fact that the E-terms embodied in
        # existing star catalogues were computed from a variety of
        # aberration constants.  Rather than adopting one of the old
        # constants the latest value is used here.
        #
        # References:
        # 1  Smith, C.A. et al., 1989.  Astr.J. 97, 265.
        # 2  Yallop, B.D. et al., 1989.  Astr.J. 97, 274.

        # Julian centuries since B1950
        T = (ep - 1950e0) * 1.00002135903e-2

        # Eccentricity
        E = 0.01673011e0 - (0.00004193e0 + 0.000000126e0 * T) * T

        # Mean obliquity
        E0 = (84404.836e0 - (46.8495e0 + (0.00319e0 + 0.00181e0 * T) * T) * T) * cls.AS2R

        # Mean longitude of perihelion
        P = (1015489.951e0 + (6190.67e0 + (1.65e0 + 0.012e0 * T) * T) * T) * cls.AS2R

        # E-terms
        EK = E * 20.49552e0 * cls.AS2R
        CP = np.cos(P)

        return np.array([EK * np.sin(P), -EK * CP * np.cos(E0), -EK * CP * np.sin(E0)])

    @classmethod
    def airmas(cls, zd):
        # +
        # - - - - - - -
        # A I R M A S
        # - - - - - - -
        #
        # Air mass at given zenith distance (double precision)
        #
        # Given:
        # ZD     d     Observed zenith distance (radians)
        #
        # The result is an estimate of the air mass, in units of that
        # at the zenith.
        #
        # Notes:
        #
        # 1)  The "observed" zenith distance referred to above means "as
        # affected by refraction".
        #
        # 2)  Uses Hardie's (1962) polynomial fit to Bemporad's data for
        # the relative air mass, X, in units of thickness at the zenith
        # as tabulated by Schoenberg (1929). This is adequate for all
        # normal needs as it is accurate to better than 0.1% up to X =
        # 6.8 and better than 1% up to X = 10. Bemporad's tabulated
        # values are unlikely to be trustworthy to such accuracy
        # because of variations in density, pressure and other
        # conditions in the atmosphere from those assumed in his work.
        #
        # 3)  The sign of the ZD is ignored.
        #
        # 4)  At zenith distances greater than about ZD = 87 degrees the
        # air mass is held constant to avoid arithmetic overflows.
        #
        # References:
        # Hardie, R.H., 1962, in "Astronomical Techniques"
        # ed. W.A. Hiltner, University of Chicago Press, p180.
        # Schoenberg, E., 1929, Hdb. d. Ap.,
        # Berlin, Julius Springer, 2, 268.

        SECZM1 = 1e0 / (np.cos(np.minimum(1.52e0, np.abs(zd)))) - 1e0
        return 1e0 + SECZM1 * (0.9981833e0 - SECZM1 * (0.002875e0 + 0.0008083e0 * SECZM1))

    @classmethod
    def altaz(cls, ha, dec, phi):

        # +
        # - - - - - -
        # A L T A Z
        # - - - - - -
        #
        # Positions, velocities and accelerations for an altazimuth
        # telescope mount.
        #
        # (double precision)
        #
        # Given:
        # HA      d     hour angle
        # DEC     d     declination
        # PHI     d     observatory latitude
        #
        # Returned:
        # A      d     azimuth
        # AD     d        "    velocity
        # ADD    d        "    acceleration
        # E      d     elevation
        # ED     d         "     velocity
        # EDD    d         "     acceleration
        # Q      d     parallactic angle
        # QD     d         "      "   velocity
        # QDD    d         "      "   acceleration
        #
        # Notes:
        #
        # 1)  Natural units are used throughout.  HA, DEC, PHI, AZ, EL
        # and ZD are in radians.  The velocities and accelerations
        # assume constant declination and constant rate of change of
        # hour angle (as for tracking a star);  the units of AZD, ELD
        # and PAD are radians per radian of HA, while the units of AZDD,
        # ELDD and PADD are radians per radian of HA squared.  To
        # convert into practical degree- and second-based units:
        #
        # angles * 360/2pi -> degrees
        # velocities * (2pi/86400)*(360/2pi) -> degree/sec
        # accelerations * ((2pi/86400)**2)*(360/2pi) -> degree/sec/sec
        #
        # Note that the seconds here are sidereal rather than SI.  One
        # sidereal second is about 0.99727 SI seconds.
        #
        # The velocity and acceleration factors assume the sidereal
        # tracking case.  Their respective numerical values are (exactly)
        # 1/240 and (approximately) 1/3300236.9.
        #
        # 2)  Azimuth is returned in the range 0-2pi;  north is zero,
        # and east is +pi/2.  Elevation and parallactic angle are
        # returned in the range +/-pi.  Parallactic angle is +ve for
        # a star west of the meridian and is the angle NP-star-zenith.
        #
        # 3)  The latitude is geodetic as opposed to geocentric.  The
        # hour angle and declination are topocentric.  Refraction and
        # deficiencies in the telescope mounting are ignored.  The
        # purpose of the routine is to give the general form of the
        # quantities.  The details of a real telescope could profoundly
        # change the results, especially close to the zenith.
        #
        # 4)  No range checking of arguments is carried out.
        #
        # 5)  In applications which involve many such calculations, rather
        # than calling the present routine it will be more efficient to
        # use inline code, having previously computed fixed terms such
        # as sine and cosine of latitude, and (for tracking a star)
        # sine and cosine of declination.

        # Useful functions
        SH = np.sin(ha)
        CH = np.cos(ha)
        SD = np.sin(dec)
        CD = np.cos(dec)
        SP = np.sin(phi)
        CP = np.cos(phi)
        CHCD = CH * CD
        SDCP = SD * CP
        X = -CHCD * SP + SDCP
        Y = -SH * CD
        Z = CHCD * CP + SD * SP
        RSQ = X * X + Y * Y
        R = np.sqrt(RSQ)

        # Azimuth and elevation
        A = np.where(RSQ == 0e0, 0e0, np.arctan2(Y, X))
        A = np.where(A < 0e0, A + cls.D2PI, A)
        E = np.arctan2(Z, R)

        # Parallactic angle
        C = CD * SP - CH * SDCP
        S = SH * CP

        Q = np.where(C * C + S * S > 0, np.arctan2(S, C), cls.DPI - ha)

        # Velocities and accelerations (clamped at zenith/nadir)
        R = np.clip(R, cls.TINY, None)
        RSQ = np.clip(RSQ, cls.TINY, None)

        QD = -X * CP / RSQ
        AD = SP + Z * QD
        ED = CP * Y / R
        EDR = ED / R
        ADD = EDR * (Z * SP + (2e0 - RSQ) * QD)
        EDD = -R * QD * AD
        QDD = EDR * (SP + 2e0 * Z * QD)

        return A, AD, ADD, E, ED, EDD, Q, QD, QDD

    @classmethod
    def amp(cls, ra, da, date, eq):
        # +
        # - - - -
        # A M P
        # - - - -
        #
        # Convert star RA,Dec from geocentric apparent to mean place
        #
        # The mean coordinate system is the post IAU 1976 system,
        # loosely called FK5.
        #
        # Given:
        # RA       d      apparent RA (radians)
        # DA       d      apparent Dec (radians)
        # DATE     d      TDB for apparent place (JD-2400000.5)
        # EQ       d      equinox:  Julian epoch of mean place
        #
        # Returned:
        # RM       d      mean RA (radians)
        # DM       d      mean Dec (radians)
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Notes:
        #
        # 1)  The distinction between the required TDB and TT is always
        # negligible.  Moreover, for all but the most critical
        # applications UTC is adequate.
        #
        # 2)  Iterative techniques are used for the aberration and light
        # deflection corrections so that the routines sla_AMP (or
        # sla_AMPQK) and sla_MAP (or sla_MAPQK) are accurate inverses;
        # even at the edge of the Sun's disc the discrepancy is only
        # about 1 nanoarcsecond.
        #
        # 3)  Where multiple apparent places are to be converted to mean
        # places, for a fixed date and equinox, it is more efficient to
        # use the sla_MAPPA routine to compute the required parameters
        # once, followed by one call to sla_AMPQK per star.
        #
        # 4)  The accuracy is sub-milliarcsecond, limited by the
        # precession-nutation model (IAU 1976 precession, Shirai &
        # Fukushima 2001 forced nutation and precession corrections).
        #
        # 5)  The accuracy is further limited by the routine sla_EVP, called
        # by sla_MAPPA, which computes the Earth position and velocity
        # using the methods of Stumpff.  The maximum error is about
        # 0.3 mas.
        #
        # Depends:  sla_MAPPA, sla_AMPQK

        amprms = cls.mappa(eq, date)
        rm, dm = cls.ampqk(ra, da, amprms)

        return rm, dm

    @classmethod
    def ampqk(cls, ra, da, amprms):
        # +
        # - - - - - -
        # A M P Q K
        # - - - - - -
        #
        # Convert star RA,Dec from geocentric apparent to mean place
        #
        # The mean coordinate system is the post IAU 1976 system,
        # loosely called FK5.
        #
        # Use of this routine is appropriate when efficiency is important
        # and where many star positions are all to be transformed for
        # one epoch and equinox.  The star-independent parameters can be
        # obtained by calling the sla_MAPPA routine.
        #
        # Given:
        # RA       d      apparent RA (radians)
        # DA       d      apparent Dec (radians)
        #
        # AMPRMS   d(21)  star-independent mean-to-apparent parameters:
        #
        # (0)      time interval for proper motion (Julian years)
        # (1)    barycentric position of the Earth (AU)
        # (2)    heliocentric direction of the Earth (unit vector)
        # (3)      (grav rad Sun)*2/(Sun-Earth distance)
        # (4)   ABV: barycentric Earth velocity in units of c
        # (5)     sqrt(1-v**2) where v=modulus(ABV)
        # (6)  precession/nutation (3,3) matrix
        #
        # Returned:
        # RM       d      mean RA (radians)
        # DM       d      mean Dec (radians)
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Note:
        #
        # Iterative techniques are used for the aberration and
        # light deflection corrections so that the routines
        # sla_AMP (or sla_AMPQK) and sla_MAP (or sla_MAPQK) are
        # accurate inverses;  even at the edge of the Sun's disc
        # the discrepancy is only about 1 nanoarcsecond.
        #
        # Depends:  sla_DCS2C, sla_DIMXV, sla_DVDV, sla_DVN, sla_DCC2S,
        # sla_DRANRM

        # Unpack scalar and vector parameters
        GR2E = amprms[3]
        AB1 = amprms[5]

        EHN = amprms[2]
        ABV = amprms[4]

        # Apparent RA,Dec to Cartesian
        P3 = cls.dcs2c(ra, da)

        # Precession and nutation
        P2 = cls.dimxv(amprms[6], P3)

        # Aberration
        AB1P1 = AB1 + 1e0
        P1 = P2.copy()

        for _ in range(2):
            P1DV = cls.dvdv(P1, ABV)
            P1DVP1 = 1e0 + P1DV
            W = 1e0 + P1DV / AB1P1
            P1 = (P1DVP1 * P2 - W * ABV) / AB1

            P3, W = cls.dvn(P1)
            P1 = P3.copy()

        # Light deflection
        P = P1.copy()

        for _ in range(5):
            PDE = cls.dvdv(P, EHN)
            PDEP1 = 1e0 + PDE
            W = PDEP1 - GR2E * PDE
            P = (PDEP1 * P1 - GR2E * EHN) / W

            P2, W = cls.dvn(P)
            P = P2.copy()

        # Mean RA,Dec
        rm, dm = cls.dcc2s(P)
        rm = cls.dranrm(rm)

        return rm, dm

    @classmethod
    def dvdv(cls, va, vb):
        # +
        # - - - - -
        # D V D V
        # - - - - -
        #
        # Scalar product of two 3-vectors
        #
        # Given:
        # VA      dp(3)     first vector
        # VB      dp(3)     second vector

        return np.dot(va, vb)

    @classmethod
    def dvn(cls, v):
        # +
        # - - - -
        # D V N
        # - - - -
        #
        # Normalizes a 3-vector also giving the modulus (double precision)
        #
        # Given:
        # V       d(3)      vector
        #
        # Returned:
        # UV      d(3)      unit vector in direction of V
        # VM      d         modulus of V
        #
        # Notes:
        # If the modulus of V is zero, UV is set to zero as well.

        # Modulus.
        vm = np.sqrt(np.dot(v, v))

        # Normalize the vector.
        vm = np.where(vm <= 0e0, 1e0, vm)

        return v / vm, vm

    @classmethod
    def dimxv(cls, dm, va):
        # +
        # - - - - - -
        # D I M X V
        # - - - - - -
        #
        # Performs the 3-D backward unitary transformation:
        #
        # vector VB = (inverse of matrix DM) * vector VA
        #
        # (double precision)
        #
        # (n.b.  the matrix must be unitary, as this routine assumes that
        # the inverse and transpose are identical)
        #
        # Given:
        # DM       dp(3,3)    matrix
        # VA       dp(3)      vector
        #
        # Returned:
        # VB       dp(3)      result vector

        # Inverse of matrix DM * vector VA -> vector VW
        return dm.I @ va

    @classmethod
    def mappa(cls, eq, date):
        # +
        # - - - - - -
        # M A P P A
        # - - - - - -
        #
        # Compute star-independent parameters in preparation for
        # conversions between mean place and geocentric apparent place.
        #
        # The parameters produced by this routine are required in the
        # parallax, light deflection, aberration, and precession/nutation
        # parts of the mean/apparent transformations.
        #
        # The reference frames and timescales used are post IAU 1976.
        #
        # Given:
        # EQ       d      epoch of mean equinox to be used (Julian)
        # DATE     d      TDB (JD-2400000.5)
        #
        # Returned:
        # AMPRMS   d(7)  star-independent mean-to-apparent parameters:
        #
        # (0)      time interval for proper motion (Julian years)
        # (1)    barycentric position of the Earth (AU)
        # (2)    heliocentric direction of the Earth (unit vector)
        # (3)      (grav rad Sun)*2/(Sun-Earth distance)
        # (4)   ABV: barycentric Earth velocity in units of c
        # (5)     sqrt(1-v**2) where v=modulus(ABV)
        # (6)  precession/nutation (3,3) matrix
        #
        # References:
        # 1984 Astronomical Almanac, pp B39-B41.
        # (also Lederle & Schwan, Astron. Astrophys. 134,
        # 1-6, 1984)
        #
        # Notes:
        #
        # 1)  For DATE, the distinction between the required TDB and TT
        # is always negligible.  Moreover, for all but the most
        # critical applications UTC is adequate.
        #
        # 2)  The vectors AMPRMS(2-4) and AMPRMS(5-7) are referred to
        # the mean equinox and equator of epoch EQ.
        #
        # 3)  The parameters AMPRMS produced by this routine are used by
        # sla_AMPQK, sla_MAPQK and sla_MAPQKZ.
        #
        # 4)  The accuracy is sub-milliarcsecond, limited by the
        # precession-nutation model (IAU 1976 precession, Shirai &
        # Fukushima 2001 forced nutation and precession corrections).
        #
        # 5)  A further limit to the accuracy of routines using the parameter
        # array AMPRMS is imposed by the routine sla_EVP, used here to
        # compute the Earth position and velocity by the methods of
        # Stumpff.  The maximum error in the resulting aberration
        # corrections is about 0.3 milliarcsecond.
        #
        # Called:
        # sla_EPJ         MDJ to Julian epoch
        # sla_EVP         earth position & velocity
        # sla_DVN         normalize vector
        # sla_PRENUT      precession/nutation matrix
        # Time interval for proper motion correction

        amprms = [None for _ in range(7)]
        amprms[0] = cls.epj(date) - eq

        # Get Earth barycentric and heliocentric position and velocity
        ebd, amprms[1], EHD, EH = cls.evp(date, eq)

        # Heliocentric direction of earth (normalized) and modulus
        amprms[2], E = cls.dvn(EH)

        # Light deflection parameter
        amprms[3] = cls.GR2 / E

        # Aberration parameters
        amprms[4] = ebd * cls.CR

        VN, VM = cls.dvn(ebd * cls.CR)
        amprms[5] = np.sqrt(1e0 - VM * VM)

        # Precession/nutation matrix
        amprms[6] = cls.prenut(eq, date)

        return amprms

    @classmethod
    def epj(cls, date):
        # +
        # - - - -
        # E P J
        # - - - -
        #
        # Conversion of Modified Julian Date to Julian Epoch (double precision)
        #
        # Given:
        # DATE     dp       Modified Julian Date (JD - 2400000.5)
        #
        # The result is the Julian Epoch.
        #
        # Reference:
        # Lieske,J.H., 1979. Astron.Astrophys.,73,282.

        return 2000e0 + (date - 51544.5e0) / 365.25

    @classmethod
    def evp(cls, date, deqx):
        # +
        # - - - -
        # E V P
        # - - - -
        #
        # Barycentric and heliocentric velocity and position of the Earth
        #
        # All arguments are double precision
        #
        # Given:
        #
        # DATE          TDB (loosely ET) as a Modified Julian Date
        # (JD-2400000.5)
        #
        # DEQX          Julian Epoch (e.g. 2000.0D0) of mean equator and
        # equinox of the vectors returned.  If DEQX .LE. 0D0,
        # all vectors are referred to the mean equator and
        # equinox (FK5) of epoch DATE.
        #
        # Returned (all 3D Cartesian vectors):
        #
        # DVB,DPB       barycentric velocity, position (AU/s, AU)
        # DVH,DPH       heliocentric velocity, position (AU/s, AU)
        #
        # Depends:  sla_EPJ, sla_PREC
        #
        # Notes:
        #
        # 1  This routine is accurate enough for many purposes but faster and
        # more compact than the sla_EPV routine.  The maximum deviations
        # from the JPL DE96 ephemeris are as follows:
        #
        # barycentric velocity         0.42  m/s
        # barycentric position         6900  km
        #
        # heliocentric velocity        0.42  m/s
        # heliocentric position        1600  km
        #
        # 2  The routine is adapted from the BARVEL and BARCOR subroutines of
        # Stumpff (1980).  Most of the changes are merely cosmetic and do
        # not affect the results at all.  However, some adjustments have
        # been made so as to give results that refer to the IAU 1976 'FK5'
        # equinox and precession, although the differences these changes
        # make relative to the results from Stumpff's original 'FK4' version
        # are smaller than the inherent accuracy of the algorithm.  One
        # minor shortcoming in the original routines that has NOT been
        # corrected is that better numerical accuracy could be achieved if
        # the various polynomial evaluations were nested.
        #
        # Reference:
        #
        # Stumpff, P., Astron.Astrophys.Suppl.Ser. 41, 1-8 (1980).

        # Constants DCFEL(I,K) of fast changing elements
        # I=1                I=2              I=3
        DCFEL = np.array([
            [1.7400353e0, 6.2833195099091e2, 5.2796e-6],
            [6.2565836e0, 6.2830194572674e2, -2.6180e-6],
            [4.7199666e0, 8.3997091449254e3, -1.9780e-5],
            [1.9636505e-1, 8.4334662911720e3, -5.6044e-5],
            [4.1547339e0, 5.2993466764997e1, 5.8845e-6],
            [4.6524223e0, 2.1354275911213e1, 5.6797e-6],
            [4.2620486e0, 7.5025342197656e0, 5.5317e-6],
            [1.4740694e0, 3.8377331909193e0, 5.6093e-6]
        ])

        #
        # Constants DCEPS and CCSEL(I,K) of slowly changing elements
        # I=1           I=2           I=3
        DCEPS = np.array([4.093198e-1, -2.271110e-4, -2.860401e-8])
        CCSEL = np.array([
            [1.675104e-2, -4.179579e-5, -1.260516e-7],
            [2.220221e-1, 2.809917e-2, 1.852532e-5],
            [1.589963, 3.418075e-2, 1.430200e-5],
            [2.994089, 2.590824e-2, 4.155840e-6],
            [8.155457e-1, 2.486352e-2, 6.836840e-6],
            [1.735614, 1.763719e-2, 6.370440e-6],
            [1.968564, 1.524020e-2, -2.517152e-6],
            [1.282417, 8.703393e-3, 2.289292e-5],
            [2.280820, 1.918010e-2, 4.484520e-6],
            [4.833473e-2, 1.641773e-4, -4.654200e-7],
            [5.589232e-2, -3.455092e-4, -7.388560e-7],
            [4.634443e-2, -2.658234e-5, 7.757000e-8],
            [8.997041e-3, 6.329728e-6, -1.939256e-9],
            [2.284178e-2, -9.941590e-5, 6.787400e-8],
            [4.350267e-2, -6.839749e-5, -2.714956e-7],
            [1.348204e-2, 1.091504e-5, 6.903760e-7],
            [3.106570e-2, -1.665665e-4, -1.590188e-7]
        ])

        #
        # Constants of the arguments of the short-period perturbations
        # by the planets:   DCARGS(I,K)
        # I=1               I=2
        DCARGS = np.array([
            [5.0974222e0, -7.8604195454652e2],
            [3.9584962e0, -5.7533848094674e2],
            [1.6338070e0, -1.1506769618935e3],
            [2.5487111e0, -3.9302097727326e2],
            [4.9255514e0, -5.8849265665348e2],
            [1.3363463e0, -5.5076098609303e2],
            [1.6072053e0, -5.2237501616674e2],
            [1.3629480e0, -1.1790629318198e3],
            [5.5657014e0, -1.0977134971135e3],
            [5.0708205e0, -1.5774000881978e2],
            [3.9318944e0, 5.2963464780000e1],
            [4.8989497e0, 3.9809289073258e1],
            [1.3097446e0, 7.7540959633708e1],
            [3.5147141e0, 7.9618578146517e1],
            [3.5413158e0, -5.4868336758022e2]
        ])

        #
        # Amplitudes CCAMPS(N,K) of the short-period perturbations
        # N=1          N=2          N=3          N=4          N=5
        CCAMPS = np.array([
            [-2.279594e-5, 1.407414e-5, 8.273188e-6, 1.340565e-5, -2.490817e-7],
            [-3.494537e-5, 2.860401e-7, 1.289448e-7, 1.627237e-5, -1.823138e-7],
            [6.593466e-7, 1.322572e-5, 9.258695e-6, -4.674248e-7, -3.646275e-7],
            [1.140767e-5, -2.049792e-5, -4.747930e-6, -2.638763e-6, -1.245408e-7],
            [9.516893e-6, -2.748894e-6, -1.319381e-6, -4.549908e-6, -1.864821e-7],
            [7.310990e-6, -1.924710e-6, -8.772849e-7, -3.334143e-6, -1.745256e-7],
            [-2.603449e-6, 7.359472e-6, 3.168357e-6, 1.119056e-6, -1.655307e-7],
            [-3.228859e-6, 1.308997e-7, 1.013137e-7, 2.403899e-6, -3.736225e-7],
            [3.442177e-7, 2.671323e-6, 1.832858e-6, -2.394688e-7, -3.478444e-7],
            [8.702406e-6, -8.421214e-6, -1.372341e-6, -1.455234e-6, -4.998479e-8],
            [-1.488378e-6, -1.251789e-5, 5.226868e-7, -2.049301e-7, 0.0],
            [-8.043059e-6, -2.991300e-6, 1.473654e-7, -3.154542e-7, 0.0],
            [3.699128e-6, -3.316126e-6, 2.901257e-7, 3.407826e-7, 0.0],
            [2.550120e-6, -1.241123e-6, 9.901116e-8, 2.210482e-7, 0.0],
            [-6.351059e-7, 2.341650e-6, 1.061492e-6, 2.878231e-7, 0.0]
        ])

        #
        # Constants of the secular perturbations in longitude
        # CCSEC3 and CCSEC(N,K)
        # N=1           N=2           N=3
        CCSEC3 = -7.757020e-8
        CCSEC = np.array([
            [1.289600e-6, 5.550147e-1, 2.076942],
            [3.102810e-5, 4.035027, 3.525565e-1],
            [9.124190e-6, 9.990265e-1, 2.622706],
            [9.793240e-7, 5.508259, 1.559103e1]
        ])

        # Sidereal rate DCSLD in longitude, rate CCSGD in mean anomaly
        DCSLD = 1.990987e-7
        CCSGD = 1.990969e-7

        # Some constants used in the calculation of the lunar contribution
        CCKM = 3.122140e-5
        CCMLD = 2.661699e-6
        CCFDI = 2.399485e-7

        #
        # Constants DCARGM(I,K) of the arguments of the perturbations
        # of the motion of the Moon
        # I=1               I=2
        DCARGM = np.array([
            [5.1679830e0, 8.3286911095275e3],
            [5.4913150e0, -7.2140632838100e3],
            [5.9598530e0, 1.5542754389685e4]
        ])

        #
        # Amplitudes CCAMPM(N,K) of the perturbations of the Moon
        # N=1          N=2           N=3           N=4
        CCAMPM = np.array([
            [1.097594e-1, 2.896773e-7, 5.450474e-2, 1.438491e-7],
            [-2.223581e-2, 5.083103e-8, 1.002548e-2, -2.291823e-8],
            [1.148966e-2, 5.658888e-8, 8.249439e-3, 4.063015e-8]
        ])

        #
        # CCPAMV(K)=A*M*DL/DT (planets), DC1MME=1-MASS(Earth+Moon)
        CCPAMV = np.array([8.326827e-11, 1.843484e-11, 1.988712e-12, 1.881276e-12])
        DC1MME = 0.99999696e0

        # CCPAM(K)=A*M(planets), CCIM=INCLINATION(Moon)
        CCPAM = np.array([4.960906e-3, 2.727436e-3, 8.392311e-4, 1.556861e-3])
        CCIM = 8.978749e-2

        #
        # EXECUTION
        # ---------

        # Control parameter IDEQ, and time arguments
        DT = (date - 15019.5e0) / 36525e0

        dt_poly = np.array([1., DT, DT * DT])
        # Values of all elements for the instant DATE
        DLOCAL = np.mod(np.dot(DCFEL, dt_poly), cls.D2PI)
        DML = DLOCAL[0]
        FORBEL = DLOCAL[1:]

        DEPS = np.mod(np.dot(DCEPS, dt_poly), cls.D2PI)
        SORBEL = np.mod(np.dot(CCSEL, dt_poly), cls.D2PI)

        # Secular perturbations in longitude
        dt_line = np.array([1., DT])
        SN = np.sin(np.mod(np.dot(CCSEC[:, 1:], dt_line), cls.D2PI))

        # Periodic perturbations of the EMB (Earth-Moon barycentre)
        PERTL = CCSEC[0, 0] * SN(1) + CCSEC[1, 0] * SN(2) + (CCSEC[2, 0] + DT * CCSEC3) * SN(3) + CCSEC[3, 0] * SN(4)

        PERTLD = 0.0
        PERTR = 0.0
        PERTRD = 0.0
        for row_ind, (args_row, amps_row) in enumerate(zip(DCARGS, CCAMPS)):
            A = np.mod(np.dot(args_row[:2], dt_line), cls.D2PI)
            COSA = np.cos(A)
            SINA = np.sin(A)
            PERTL += amps_row[0] * COSA + amps_row[1] * SINA
            PERTR += amps_row[2] * COSA + amps_row[3] * SINA
            if row_ind < 11:
                PERTLD += (amps_row[1] * COSA - amps_row[0] * SINA) * amps_row[4]
                PERTRD += (amps_row[3] * COSA - amps_row[2] * SINA) * amps_row[4]

        # Elliptic part of the motion of the EMB
        E = SORBEL[0]
        G = FORBEL[0]
        ESQ = E * E
        DPARAM = 1e0 - ESQ
        PARAM = DPARAM
        TWOE = E + E
        TWOG = G + G
        PHI = TWOE * ((1.0 - ESQ * 0.125) * np.sin(G) + E * 0.625 * np.sin(TWOG) + ESQ * 0.54166667 * np.sin(G + TWOG))

        F = G + PHI
        SINF = np.sin(F)
        COSF = np.cos(F)
        DPSI = DPARAM / (1e0 + (E * COSF))
        PHID = TWOE * CCSGD * ((1.0 + ESQ * 1.5) * COSF + E * (1.25 - SINF * SINF * 0.5))
        PSID = CCSGD * E * SINF / np.sqrt(PARAM)

        # Perturbed heliocentric motion of the EMB
        D1PDRO = 1e0 + PERTR
        DRD = D1PDRO * (PSID + DPSI * PERTRD)
        DRLD = D1PDRO * DPSI * (DCSLD + PHID + PERTLD)
        DTL = np.mod(DML + PHI + PERTL, cls.D2PI)
        DSINLS = np.sin(DTL)
        DCOSLS = np.cos(DTL)
        DXHD = DRD * DCOSLS - DRLD * DSINLS
        DYHD = DRD * DSINLS + DRLD * DCOSLS

        # Influence of eccentricity, evection and variation on the
        # geocentric motion of the Moon
        PERTL = 0.0
        PERTLD = 0.0
        PERTP = 0.0
        PERTPD = 0.0
        for arg_row, amp_row in zip(DCARGM, CCAMPM):
            A = np.mod(np.dot(arg_row, dt_line), cls.D2PI)
            SINA = np.sin(A)
            COSA = np.cos(A)
            PERTL += amp_row[0] * SINA
            PERTLD += amp_row[1] * COSA
            PERTP += amp_row[2] * COSA
            PERTPD -= amp_row[3] * SINA

        # Heliocentric motion of the Earth
        TL = FORBEL[1] + PERTL
        SINLM = np.sin(TL)
        COSLM = np.cos(TL)
        SIGMA = CCKM / (1.0 + PERTP)
        A = SIGMA * (CCMLD + PERTLD)
        B = SIGMA * PERTPD
        DXHD = DXHD + (A * SINLM) + (B * COSLM)
        DYHD = DYHD - (A * COSLM) + (B * SINLM)
        DZHD = -(SIGMA * CCFDI * np.cos(FORBEL[2]))

        # Barycentric motion of the Earth
        DXBD = DXHD * DC1MME
        DYBD = DYHD * DC1MME
        DZBD = DZHD * DC1MME

        SINLP = np.zeros(4)
        COSLP = np.zeros(4)
        for K in range(4):
            PLON = FORBEL[K + 3]
            POMG = SORBEL[K + 1]
            PECC = SORBEL[K + 9]
            TL = np.mod(PLON + 2.0 * PECC * np.sin(PLON - POMG), cls.D2PI)
            SINLP[K] = np.sin(TL)
            COSLP[K] = np.cos(TL)
            DXBD = DXBD + (CCPAMV[K] * (SINLP[K] + PECC * np.sin(POMG)))
            DYBD = DYBD - (CCPAMV[K] * (COSLP[K] + PECC * np.cos(POMG)))
            DZBD = DZBD - (CCPAMV[K] * SORBEL[K + 13] * np.cos(PLON - SORBEL[K + 5]))

        # Transition to mean equator of date
        DCOSEP = np.cos(DEPS)
        DSINEP = np.sin(DEPS)
        DYAHD = DCOSEP * DYHD - DSINEP * DZHD
        DZAHD = DSINEP * DYHD + DCOSEP * DZHD
        DYABD = DCOSEP * DYBD - DSINEP * DZBD
        DZABD = DSINEP * DYBD + DCOSEP * DZBD

        # Heliocentric coordinates of the Earth
        DR = DPSI * D1PDRO
        FLATM = CCIM * np.sin(FORBEL(3))
        A = SIGMA * np.cos(FLATM)
        DXH = DR * DCOSLS - (A * COSLM)
        DYH = DR * DSINLS - (A * SINLM)
        DZH = -(SIGMA * np.sin(FLATM))

        # Barycentric coordinates of the Earth
        DXB = DXH * DC1MME
        DYB = DYH * DC1MME
        DZB = DZH * DC1MME
        for K in range(4):
            FLAT = SORBEL[K + 13] * np.sin(FORBEL[K + 3] - SORBEL[K + 5])
            A = CCPAM[K] * (1.0 - SORBEL[K + 9] * np.cos(FORBEL[K + 3] - SORBEL[K + 1]))
            B = A * np.cos(FLAT)
            DXB = DXB - (B * COSLP[K])
            DYB = DYB - (B * SINLP[K])
            DZB = DZB - (A * np.sin(FLAT))

        # Transition to mean equator of date
        DYAH = DCOSEP * DYH - DSINEP * DZH
        DZAH = DSINEP * DYH + DCOSEP * DZH
        DYAB = DCOSEP * DYB - DSINEP * DZB
        DZAB = DSINEP * DYB + DCOSEP * DZB

        # Copy result components into vectors, correcting for FK4 equinox
        DEPJ = cls.epj(date)
        DEQCOR = cls.AS2R * (0.035e0 + 0.00085e0 * (DEPJ - cls.B1950))
        dvh = np.array([DXHD - DEQCOR * DYAHD, DYAHD + DEQCOR * DXHD, DZAHD])
        dvb = np.array([DXBD - DEQCOR * DYABD, DYABD + DEQCOR * DXBD, DZABD])
        dph = np.array([DXH - DEQCOR * DYAH, DYAH + DEQCOR * DXH, DZAH])
        dpb = np.array([DXB - DEQCOR * DYAB, DYAB + DEQCOR * DXB, DZAB])

        # Was precession to another equinox requested?
        if deqx > 0.:
            # Yes: compute precession matrix from MJD DATE to Julian epoch DEQX
            DPREMA = cls.prec(DEPJ, deqx)

            # Rotate DVH
            dvh = DPREMA @ dvh

            # Rotate DVB
            dvb = DPREMA @ dvb

            # Rotate DPH
            dph = DPREMA @ dph

            # Rotate DPB
            dpb = DPREMA @ dpb

        return dvb, dpb, dvh, dph

    @classmethod
    def prenut(cls, epoch, date):
        # +
        # - - - - - - -
        # P R E N U T
        # - - - - - - -
        #
        # Form the matrix of precession and nutation (SF2001)
        # (double precision)
        #
        # Given:
        # EPOCH   dp         Julian Epoch for mean coordinates
        # DATE    dp         Modified Julian Date (JD-2400000.5)
        # for true coordinates
        #
        # Returned:
        # RMATPN  dp(3,3)    combined precession/nutation matrix
        #
        # Depends:  sla_PREC, sla_EPJ, sla_NUT, sla_DMXM
        #
        # Notes:
        #
        # 1)  The epoch and date are TDB (loosely ET).  TT will do, or even
        # UTC.
        #
        # 2)  The matrix is in the sense   V(true) = RMATPN * V(mean)
        # -

        # Precession
        RMATP = cls.prec(epoch, cls.epj(date))

        # Nutation
        RMATN = cls.nut(date)

        return cls.dmxm(RMATN, RMATP)

    @classmethod
    def dmxm(cls, a, b):
        # +
        # - - - - -
        # D M X M
        # - - - - -
        #
        # Product of two 3x3 matrices:
        #
        # matrix C  =  matrix A  x  matrix B
        #
        # (double precision)
        #
        # Given:
        # A      dp(3,3)        matrix
        # B      dp(3,3)        matrix
        #
        # Returned:
        # C      dp(3,3)        matrix result
        # -

        # Multiply into scratch matrix
        return a @ b

    @classmethod
    def nut(cls, date):
        # +
        # - - - -
        # N U T
        # - - - -
        #
        # Form the matrix of nutation for a given date - Shirai & Fukushima
        # 2001 theory (double precision)
        #
        # Reference:
        # Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
        #
        # Given:
        # DATE    d          TDB (loosely ET) as Modified Julian Date
        # (=JD-2400000.5)
        # Returned:
        # RMATN   d(3,3)     nutation matrix
        #
        # Notes:
        #
        # 1  The matrix is in the sense  v(true) = rmatn * v(mean) .
        # where v(true) is the star vector relative to the true equator and
        # equinox of date and v(mean) is the star vector relative to the
        # mean equator and equinox of date.
        #
        # 2  The matrix represents forced nutation (but not free core
        # nutation) plus corrections to the IAU~1976 precession model.
        #
        # 3  Earth attitude predictions made by combining the present nutation
        # matrix with IAU~1976 precession are accurate to 1~mas (with
        # respect to the ICRS) for a few decades around 2000.
        #
        # 4  The distinction between the required TDB and TT is always
        # negligible.  Moreover, for all but the most critical applications
        # UTC is adequate.
        #
        # Depends:   sla_NUTC, sla_DEULER
        # -

        # Nutation components and mean obliquity
        DPSI, DEPS, EPS0 = cls.nutc(date)

        return cls.deuler('XZX', EPS0, -DPSI, -(EPS0 + DEPS))

    @classmethod
    def nutc(cls, date):
        # +
        # - - - - -
        # N U T C
        # - - - - -
        #
        # Nutation:  longitude & obliquity components and mean obliquity,
        # using the Shirai & Fukushima (2001) theory.
        #
        # Given:
        # DATE        d    TDB (loosely ET) as Modified Julian Date
        # (JD-2400000.5)
        # Returned:
        # DPSI,DEPS   d    nutation in longitude,obliquity
        # EPS0        d    mean obliquity
        #
        # Notes:
        #
        # 1  The routine predicts forced nutation (but not free core nutation)
        # plus corrections to the IAU 1976 precession model.
        #
        # 2  Earth attitude predictions made by combining the present nutation
        # model with IAU 1976 precession are accurate to 1 mas (with respect
        # to the ICRF) for a few decades around 2000.
        #
        # 3  The sla_NUTC80 routine is the equivalent of the present routine
        # but using the IAU 1980 nutation theory.  The older theory is less
        # accurate, leading to errors as large as 350 mas over the interval
        # 1900-2100, mainly because of the error in the IAU 1976 precession.
        #
        # References:
        #
        # Shirai, T. & Fukushima, T., Astron.J. 121, 3270-3283 (2001).
        #
        # Fukushima, T., Astron.Astrophys. 244, L11 (1991).
        #
        # Simon, J. L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
        # Francou, G. & Laskar, J., Astron.Astrophys. 282, 663 (1994).
        # -

        # The SF2001 forced nutation model

        # Coefficients of fundamental angles
        NA = np.array([
            [0, 0, 0, 0, -1, 0, 0, 0, 0],
            [0, 0, 2, -2, 2, 0, 0, 0, 0],
            [0, 0, 2, 0, 2, 0, 0, 0, 0],
            [0, 0, 0, 0, -2, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 2, -2, 2, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 0, 1, 0, 0, 0, 0],
            [1, 0, 2, 0, 2, 0, 0, 0, 0],
            [0, -1, 2, -2, 2, 0, 0, 0, 0],
            [0, 0, 2, -2, 1, 0, 0, 0, 0],
            [-1, 0, 2, 0, 2, 0, 0, 0, 0],
            [-1, 0, 0, 2, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, -1, 0, 0, 0, 0],
            [-1, 0, 2, 2, 2, 0, 0, 0, 0],
            [1, 0, 2, 0, 1, 0, 0, 0, 0],
            [-2, 0, 2, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 2, 0, 0, 0, 0, 0],
            [0, 0, 2, 2, 2, 0, 0, 0, 0],
            [2, 0, 0, -2, 0, 0, 0, 0, 0],
            [2, 0, 2, 0, 2, 0, 0, 0, 0],
            [1, 0, 2, -2, 2, 0, 0, 0, 0],
            [-1, 0, 2, 0, 1, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 1, 0, 0, 0, 0],
            [-1, 0, 0, 2, 1, 0, 0, 0, 0],
            [0, 2, 2, -2, 2, 0, 0, 0, 0],
            [0, 0, 2, -2, 0, 0, 0, 0, 0],
            [-1, 0, 0, 2, -1, 0, 0, 0, 0],
            [0, 1, 0, 0, -1, 0, 0, 0, 0],
            [0, 2, 0, 0, 0, 0, 0, 0, 0],
            [-1, 0, 2, 2, 1, 0, 0, 0, 0],
            [1, 0, 2, 2, 2, 0, 0, 0, 0],
            [0, 1, 2, 0, 2, 0, 0, 0, 0],
            [-2, 0, 2, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 2, 1, 0, 0, 0, 0],
            [0, -1, 2, 0, 2, 0, 0, 0, 0],
            [0, 0, 0, 2, 1, 0, 0, 0, 0],
            [1, 0, 2, -2, 1, 0, 0, 0, 0],
            [2, 0, 0, -2, -1, 0, 0, 0, 0],
            [2, 0, 2, -2, 2, 0, 0, 0, 0],
            [2, 0, 2, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 2, -1, 0, 0, 0, 0],
            [0, -1, 2, -2, 1, 0, 0, 0, 0],
            [-1, -1, 0, 2, 0, 0, 0, 0, 0],
            [2, 0, 0, -2, 1, 0, 0, 0, 0],
            [1, 0, 0, 2, 0, 0, 0, 0, 0],
            [0, 1, 2, -2, 1, 0, 0, 0, 0],
            [1, -1, 0, 0, 0, 0, 0, 0, 0],
            [-2, 0, 2, 0, 2, 0, 0, 0, 0],
            [0, -1, 0, 2, 0, 0, 0, 0, 0],
            [3, 0, 2, 0, 2, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0],
            [1, -1, 2, 0, 2, 0, 0, 0, 0],
            [1, 0, 0, -1, 0, 0, 0, 0, 0],
            [-1, -1, 2, 2, 2, 0, 0, 0, 0],
            [-1, 0, 2, 0, 0, 0, 0, 0, 0],
            [2, 0, 0, 0, -1, 0, 0, 0, 0],
            [0, -1, 2, 2, 2, 0, 0, 0, 0],
            [1, 1, 2, 0, 2, 0, 0, 0, 0],
            [2, 0, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, -2, 2, -1, 0, 0, 0, 0],
            [1, 0, 2, 0, 0, 0, 0, 0, 0],
            [-1, 1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 2, 0, 0, 0, 0],
            [-1, 0, 1, 0, 1, 0, 0, 0, 0],
            [0, 0, 2, 1, 2, 0, 0, 0, 0],
            [-1, 1, 0, 1, 1, 0, 0, 0, 0],
            [-1, 0, 2, 4, 2, 0, 0, 0, 0],
            [0, -2, 2, -2, 1, 0, 0, 0, 0],
            [1, 0, 2, 2, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, -2, 0, 0, 0, 0],
            [-2, 0, 2, 2, 2, 0, 0, 0, 0],
            [1, 1, 2, -2, 2, 0, 0, 0, 0],
            [-2, 0, 2, 4, 2, 0, 0, 0, 0],
            [-1, 0, 4, 0, 2, 0, 0, 0, 0],
            [2, 0, 2, -2, 1, 0, 0, 0, 0],
            [1, 0, 0, -1, -1, 0, 0, 0, 0],
            [2, 0, 2, 2, 2, 0, 0, 0, 0],
            [1, 0, 0, 2, 1, 0, 0, 0, 0],
            [3, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, -2, -1, 0, 0, 0, 0],
            [3, 0, 2, -2, 2, 0, 0, 0, 0],
            [0, 0, 4, -2, 2, 0, 0, 0, 0],
            [-1, 0, 0, 4, 0, 0, 0, 0, 0],
            [0, 1, 2, 0, 1, 0, 0, 0, 0],
            [0, 0, 2, -2, 3, 0, 0, 0, 0],
            [-2, 0, 0, 4, 0, 0, 0, 0, 0],
            [-1, -1, 0, 2, 1, 0, 0, 0, 0],
            [-2, 0, 2, 0, -1, 0, 0, 0, 0],
            [0, 0, 2, 0, -1, 0, 0, 0, 0],
            [0, -1, 2, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 2, 0, 0, 0, 0],
            [0, 0, 2, -1, 2, 0, 0, 0, 0],
            [2, 1, 0, -2, 0, 0, 0, 0, 0],
            [0, 0, 2, 4, 2, 0, 0, 0, 0],
            [-1, -1, 0, 2, -1, 0, 0, 0, 0],
            [-1, 1, 0, 2, 0, 0, 0, 0, 0],
            [1, -1, 0, 0, 1, 0, 0, 0, 0],
            [0, -1, 2, -2, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, -2, 0, 0, 0, 0],
            [1, -1, 2, 2, 2, 0, 0, 0, 0],
            [1, 0, 0, 2, -1, 0, 0, 0, 0],
            [-1, 1, 2, 2, 2, 0, 0, 0, 0],
            [3, 0, 2, 0, 1, 0, 0, 0, 0],
            [0, 1, 2, 2, 2, 0, 0, 0, 0],
            [1, 0, 2, -2, 0, 0, 0, 0, 0],
            [-1, 0, -2, 4, -1, 0, 0, 0, 0],
            [-1, -1, 2, 2, 1, 0, 0, 0, 0],
            [0, -1, 2, 2, 1, 0, 0, 0, 0],
            [2, -1, 2, 0, 2, 0, 0, 0, 0],
            [0, 0, 0, 2, 2, 0, 0, 0, 0],
            [1, -1, 2, 0, 1, 0, 0, 0, 0],
            [-1, 1, 2, 0, 2, 0, 0, 0, 0],
            [0, 1, 0, 2, 0, 0, 0, 0, 0],
            [0, 1, 2, -2, 0, 0, 0, 0, 0],
            [0, 3, 2, -2, 2, 0, 0, 0, 0],
            [0, 0, 0, 1, 1, 0, 0, 0, 0],
            [-1, 0, 2, 2, 0, 0, 0, 0, 0],
            [2, 1, 2, 0, 2, 0, 0, 0, 0],
            [1, 1, 0, 0, 1, 0, 0, 0, 0],
            [2, 0, 0, 2, 0, 0, 0, 0, 0],
            [1, 1, 2, 0, 1, 0, 0, 0, 0],
            [-1, 0, 0, 2, 2, 0, 0, 0, 0],
            [1, 0, -2, 2, 0, 0, 0, 0, 0],
            [0, -1, 0, 2, -1, 0, 0, 0, 0],
            [-1, 0, 1, 0, 2, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, -2, 2, -2, 0, 0, 0, 0],
            [0, 0, 0, 1, -1, 0, 0, 0, 0],
            [1, -1, 0, 0, -1, 0, 0, 0, 0],
            [0, 0, 0, 4, 0, 0, 0, 0, 0],
            [1, -1, 0, 2, 0, 0, 0, 0, 0],
            [1, 0, 2, 1, 2, 0, 0, 0, 0],
            [1, 0, 2, -1, 2, 0, 0, 0, 0],
            [-1, 0, 0, 2, -2, 0, 0, 0, 0],
            [0, 0, 2, 1, 1, 0, 0, 0, 0],
            [-1, 0, 2, 0, -1, 0, 0, 0, 0],
            [-1, 0, 2, 4, 1, 0, 0, 0, 0],
            [0, 0, 2, 2, 0, 0, 0, 0, 0],
            [1, 1, 2, -2, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 1, 0, 0, 0, 0],
            [-1, 0, 2, -1, 1, 0, 0, 0, 0],
            [-2, 0, 2, 2, 1, 0, 0, 0, 0],
            [2, -1, 0, 0, 0, 0, 0, 0, 0],
            [4, 0, 2, 0, 2, 0, 0, 0, 0],
            [2, 1, 2, -2, 2, 0, 0, 0, 0],
            [0, 1, 2, 1, 2, 0, 0, 0, 0],
            [1, 0, 4, -2, 2, 0, 0, 0, 0],
            [1, 1, 0, 0, -1, 0, 0, 0, 0],
            [-2, 0, 2, 4, 1, 0, 0, 0, 0],
            [2, 0, 2, 0, 0, 0, 0, 0, 0],
            [-1, 0, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 1, 0, 2, 1, 0, 0, 0, 0],
            [-1, 0, 4, 0, 1, 0, 0, 0, 0],
            [-1, 0, 0, 4, 1, 0, 0, 0, 0],
            [2, 0, 2, 2, 1, 0, 0, 0, 0],
            [2, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 5, -5, 5, -3, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 2, 0],
            [0, 0, 1, -1, 1, 0, 0, -1, 0],
            [0, 0, -1, 1, -1, 1, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 2, 0, 0],
            [0, 0, 3, -3, 3, 0, 0, -1, 0],
            [0, 0, -8, 8, -7, 5, 0, 0, 0],
            [0, 0, -1, 1, -1, 0, 2, 0, 0],
            [0, 0, -2, 2, -2, 2, 0, 0, 0],
            [0, 0, -6, 6, -6, 4, 0, 0, 0],
            [0, 0, -2, 2, -2, 0, 8, -3, 0],
            [0, 0, 6, -6, 6, 0, -8, 3, 0],
            [0, 0, 4, -4, 4, -2, 0, 0, 0],
            [0, 0, -3, 3, -3, 2, 0, 0, 0],
            [0, 0, 4, -4, 3, 0, -8, 3, 0],
            [0, 0, -4, 4, -5, 0, 8, -3, 0],
            [0, 0, 0, 0, 0, 2, 0, 0, 0],
            [0, 0, -4, 4, -4, 3, 0, 0, 0],
            [0, 1, -1, 1, -1, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 1, -1, 1, 1, 0, 0, 0],
            [0, 0, 2, -2, 2, 0, -2, 0, 0],
            [0, -1, -7, 7, -7, 5, 0, 0, 0],
            [-2, 0, 2, 0, 2, 0, 0, -2, 0],
            [-2, 0, 2, 0, 1, 0, 0, -3, 0],
            [0, 0, 2, -2, 2, 0, 0, -2, 0],
            [0, 0, 1, -1, 1, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 2],
            [0, 0, 0, 0, 0, 0, 0, 0, 1],
            [2, 0, -2, 0, -2, 0, 0, 3, 0],
            [0, 0, 1, -1, 1, 0, 0, -2, 0],
            [0, 0, -7, 7, -7, 5, 0, 0, 0]
        ])

        # Nutation series: longitude
        PSI = np.array([
            [3341.5e0, 17206241.8e0, 3.1e0, 17409.5e0],
            [-1716.8e0, -1317185.3e0, 1.4e0, -156.8e0],
            [285.7e0, -227667.0e0, 0.3e0, -23.5e0],
            [-68.6e0, -207448.0e0, 0.0e0, -21.4e0],
            [950.3e0, 147607.9e0, -2.3e0, -355.0e0],
            [-66.7e0, -51689.1e0, 0.2e0, 122.6e0],
            [-108.6e0, 71117.6e0, 0.0e0, 7.0e0],
            [35.6e0, -38740.2e0, 0.1e0, -36.2e0],
            [85.4e0, -30127.6e0, 0.0e0, -3.1e0],
            [9.0e0, 21583.0e0, 0.1e0, -50.3e0],
            [22.1e0, 12822.8e0, 0.0e0, 13.3e0],
            [3.4e0, 12350.8e0, 0.0e0, 1.3e0],
            [-21.1e0, 15699.4e0, 0.0e0, 1.6e0],
            [4.2e0, 6313.8e0, 0.0e0, 6.2e0],
            [-22.8e0, 5796.9e0, 0.0e0, 6.1e0],
            [15.7e0, -5961.1e0, 0.0e0, -0.6e0],
            [13.1e0, -5159.1e0, 0.0e0, -4.6e0],
            [1.8e0, 4592.7e0, 0.0e0, 4.5e0],
            [-17.5e0, 6336.0e0, 0.0e0, 0.7e0],
            [16.3e0, -3851.1e0, 0.0e0, -0.4e0],
            [-2.8e0, 4771.7e0, 0.0e0, 0.5e0],
            [13.8e0, -3099.3e0, 0.0e0, -0.3e0],
            [0.2e0, 2860.3e0, 0.0e0, 0.3e0],
            [1.4e0, 2045.3e0, 0.0e0, 2.0e0],
            [-8.6e0, 2922.6e0, 0.0e0, 0.3e0],
            [-7.7e0, 2587.9e0, 0.0e0, 0.2e0],
            [8.8e0, -1408.1e0, 0.0e0, 3.7e0],
            [1.4e0, 1517.5e0, 0.0e0, 1.5e0],
            [-1.9e0, -1579.7e0, 0.0e0, 7.7e0],
            [1.3e0, -2178.6e0, 0.0e0, -0.2e0],
            [-4.8e0, 1286.8e0, 0.0e0, 1.3e0],
            [6.3e0, 1267.2e0, 0.0e0, -4.0e0],
            [-1.0e0, 1669.3e0, 0.0e0, -8.3e0],
            [2.4e0, -1020.0e0, 0.0e0, -0.9e0],
            [4.5e0, -766.9e0, 0.0e0, 0.0e0],
            [-1.1e0, 756.5e0, 0.0e0, -1.7e0],
            [-1.4e0, -1097.3e0, 0.0e0, -0.5e0],
            [2.6e0, -663.0e0, 0.0e0, -0.6e0],
            [0.8e0, -714.1e0, 0.0e0, 1.6e0],
            [0.4e0, -629.9e0, 0.0e0, -0.6e0],
            [0.3e0, 580.4e0, 0.0e0, 0.6e0],
            [-1.6e0, 577.3e0, 0.0e0, 0.5e0],
            [-0.9e0, 644.4e0, 0.0e0, 0.0e0],
            [2.2e0, -534.0e0, 0.0e0, -0.5e0],
            [-2.5e0, 493.3e0, 0.0e0, 0.5e0],
            [-0.1e0, -477.3e0, 0.0e0, -2.4e0],
            [-0.9e0, 735.0e0, 0.0e0, -1.7e0],
            [0.7e0, 406.2e0, 0.0e0, 0.4e0],
            [-2.8e0, 656.9e0, 0.0e0, 0.0e0],
            [0.6e0, 358.0e0, 0.0e0, 2.0e0],
            [-0.7e0, 472.5e0, 0.0e0, -1.1e0],
            [-0.1e0, -300.5e0, 0.0e0, 0.0e0],
            [-1.2e0, 435.1e0, 0.0e0, -1.0e0],
            [1.8e0, -289.4e0, 0.0e0, 0.0e0],
            [0.6e0, -422.6e0, 0.0e0, 0.0e0],
            [0.8e0, -287.6e0, 0.0e0, 0.6e0],
            [-38.6e0, -392.3e0, 0.0e0, 0.0e0],
            [0.7e0, -281.8e0, 0.0e0, 0.6e0],
            [0.6e0, -405.7e0, 0.0e0, 0.0e0],
            [-1.2e0, 229.0e0, 0.0e0, 0.2e0],
            [1.1e0, -264.3e0, 0.0e0, 0.5e0],
            [-0.7e0, 247.9e0, 0.0e0, -0.5e0],
            [-0.2e0, 218.0e0, 0.0e0, 0.2e0],
            [0.6e0, -339.0e0, 0.0e0, 0.8e0],
            [-0.7e0, 198.7e0, 0.0e0, 0.2e0],
            [-1.5e0, 334.0e0, 0.0e0, 0.0e0],
            [0.1e0, 334.0e0, 0.0e0, 0.0e0],
            [-0.1e0, -198.1e0, 0.0e0, 0.0e0],
            [-106.6e0, 0.0e0, 0.0e0, 0.0e0],
            [-0.5e0, 165.8e0, 0.0e0, 0.0e0],
            [0.0e0, 134.8e0, 0.0e0, 0.0e0],
            [0.9e0, -151.6e0, 0.0e0, 0.0e0],
            [0.0e0, -129.7e0, 0.0e0, 0.0e0],
            [0.8e0, -132.8e0, 0.0e0, -0.1e0],
            [0.5e0, -140.7e0, 0.0e0, 0.0e0],
            [-0.1e0, 138.4e0, 0.0e0, 0.0e0],
            [0.0e0, 129.0e0, 0.0e0, -0.3e0],
            [0.5e0, -121.2e0, 0.0e0, 0.0e0],
            [-0.3e0, 114.5e0, 0.0e0, 0.0e0],
            [-0.1e0, 101.8e0, 0.0e0, 0.0e0],
            [-3.6e0, -101.9e0, 0.0e0, 0.0e0],
            [0.8e0, -109.4e0, 0.0e0, 0.0e0],
            [0.2e0, -97.0e0, 0.0e0, 0.0e0],
            [-0.7e0, 157.3e0, 0.0e0, 0.0e0],
            [0.2e0, -83.3e0, 0.0e0, 0.0e0],
            [-0.3e0, 93.3e0, 0.0e0, 0.0e0],
            [-0.1e0, 92.1e0, 0.0e0, 0.0e0],
            [-0.5e0, 133.6e0, 0.0e0, 0.0e0],
            [-0.1e0, 81.5e0, 0.0e0, 0.0e0],
            [0.0e0, 123.9e0, 0.0e0, 0.0e0],
            [-0.3e0, 128.1e0, 0.0e0, 0.0e0],
            [0.1e0, 74.1e0, 0.0e0, -0.3e0],
            [-0.2e0, -70.3e0, 0.0e0, 0.0e0],
            [-0.4e0, 66.6e0, 0.0e0, 0.0e0],
            [0.1e0, -66.7e0, 0.0e0, 0.0e0],
            [-0.7e0, 69.3e0, 0.0e0, -0.3e0],
            [0.0e0, -70.4e0, 0.0e0, 0.0e0],
            [-0.1e0, 101.5e0, 0.0e0, 0.0e0],
            [0.5e0, -69.1e0, 0.0e0, 0.0e0],
            [-0.2e0, 58.5e0, 0.0e0, 0.2e0],
            [0.1e0, -94.9e0, 0.0e0, 0.2e0],
            [0.0e0, 52.9e0, 0.0e0, -0.2e0],
            [0.1e0, 86.7e0, 0.0e0, -0.2e0],
            [-0.1e0, -59.2e0, 0.0e0, 0.2e0],
            [0.3e0, -58.8e0, 0.0e0, 0.1e0],
            [-0.3e0, 49.0e0, 0.0e0, 0.0e0],
            [-0.2e0, 56.9e0, 0.0e0, -0.1e0],
            [0.3e0, -50.2e0, 0.0e0, 0.0e0],
            [-0.2e0, 53.4e0, 0.0e0, -0.1e0],
            [0.1e0, -76.5e0, 0.0e0, 0.0e0],
            [-0.2e0, 45.3e0, 0.0e0, 0.0e0],
            [0.1e0, -46.8e0, 0.0e0, 0.0e0],
            [0.2e0, -44.6e0, 0.0e0, 0.0e0],
            [0.2e0, -48.7e0, 0.0e0, 0.0e0],
            [0.1e0, -46.8e0, 0.0e0, 0.0e0],
            [0.1e0, -42.0e0, 0.0e0, 0.0e0],
            [0.0e0, 46.4e0, 0.0e0, -0.1e0],
            [0.2e0, -67.3e0, 0.0e0, 0.1e0],
            [0.0e0, -65.8e0, 0.0e0, 0.2e0],
            [-0.1e0, -43.9e0, 0.0e0, 0.3e0],
            [0.0e0, -38.9e0, 0.0e0, 0.0e0],
            [-0.3e0, 63.9e0, 0.0e0, 0.0e0],
            [-0.2e0, 41.2e0, 0.0e0, 0.0e0],
            [0.0e0, -36.1e0, 0.0e0, 0.2e0],
            [-0.3e0, 58.5e0, 0.0e0, 0.0e0],
            [-0.1e0, 36.1e0, 0.0e0, 0.0e0],
            [0.0e0, -39.7e0, 0.0e0, 0.0e0],
            [0.1e0, -57.7e0, 0.0e0, 0.0e0],
            [-0.2e0, 33.4e0, 0.0e0, 0.0e0],
            [36.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-0.1e0, 55.7e0, 0.0e0, -0.1e0],
            [0.1e0, -35.4e0, 0.0e0, 0.0e0],
            [0.1e0, -31.0e0, 0.0e0, 0.0e0],
            [-0.1e0, 30.1e0, 0.0e0, 0.0e0],
            [-0.3e0, 49.2e0, 0.0e0, 0.0e0],
            [-0.2e0, 49.1e0, 0.0e0, 0.0e0],
            [-0.1e0, 33.6e0, 0.0e0, 0.0e0],
            [0.1e0, -33.5e0, 0.0e0, 0.0e0],
            [0.1e0, -31.0e0, 0.0e0, 0.0e0],
            [-0.1e0, 28.0e0, 0.0e0, 0.0e0],
            [0.1e0, -25.2e0, 0.0e0, 0.0e0],
            [0.1e0, -26.2e0, 0.0e0, 0.0e0],
            [-0.2e0, 41.5e0, 0.0e0, 0.0e0],
            [0.0e0, 24.5e0, 0.0e0, 0.1e0],
            [-16.2e0, 0.0e0, 0.0e0, 0.0e0],
            [0.0e0, -22.3e0, 0.0e0, 0.0e0],
            [0.0e0, 23.1e0, 0.0e0, 0.0e0],
            [-0.1e0, 37.5e0, 0.0e0, 0.0e0],
            [0.2e0, -25.7e0, 0.0e0, 0.0e0],
            [0.0e0, 25.2e0, 0.0e0, 0.0e0],
            [0.1e0, -24.5e0, 0.0e0, 0.0e0],
            [-0.1e0, 24.3e0, 0.0e0, 0.0e0],
            [0.1e0, -20.7e0, 0.0e0, 0.0e0],
            [0.1e0, -20.8e0, 0.0e0, 0.0e0],
            [-0.2e0, 33.4e0, 0.0e0, 0.0e0],
            [32.9e0, 0.0e0, 0.0e0, 0.0e0],
            [0.1e0, -32.6e0, 0.0e0, 0.0e0],
            [0.0e0, 19.9e0, 0.0e0, 0.0e0],
            [-0.1e0, 19.6e0, 0.0e0, 0.0e0],
            [0.0e0, -18.7e0, 0.0e0, 0.0e0],
            [0.1e0, -19.0e0, 0.0e0, 0.0e0],
            [0.1e0, -28.6e0, 0.0e0, 0.0e0],
            [4.0e0, 178.8e0, -11.8e0, 0.3e0],
            [39.8e0, -107.3e0, -5.6e0, -1.0e0],
            [9.9e0, 164.0e0, -4.1e0, 0.1e0],
            [-4.8e0, -135.3e0, -3.4e0, -0.1e0],
            [50.5e0, 75.0e0, 1.4e0, -1.2e0],
            [-1.1e0, -53.5e0, 1.3e0, 0.0e0],
            [-45.0e0, -2.4e0, -0.4e0, 6.6e0],
            [-11.5e0, -61.0e0, -0.9e0, 0.4e0],
            [4.4e0, -68.4e0, -3.4e0, 0.0e0],
            [7.7e0, -47.1e0, -4.7e0, -1.0e0],
            [-42.9e0, -12.6e0, -1.2e0, 4.2e0],
            [-42.8e0, 12.7e0, -1.2e0, -4.2e0],
            [-7.6e0, -44.1e0, 2.1e0, -0.5e0],
            [-64.1e0, 1.7e0, 0.2e0, 4.5e0],
            [36.4e0, -10.4e0, 1.0e0, 3.5e0],
            [35.6e0, 10.2e0, 1.0e0, -3.5e0],
            [-1.7e0, 39.5e0, 2.0e0, 0.0e0],
            [50.9e0, -8.2e0, -0.8e0, -5.0e0],
            [0.0e0, 52.3e0, 1.2e0, 0.0e0],
            [-42.9e0, -17.8e0, 0.4e0, 0.0e0],
            [2.6e0, 34.3e0, 0.8e0, 0.0e0],
            [-0.8e0, -48.6e0, 2.4e0, -0.1e0],
            [-4.9e0, 30.5e0, 3.7e0, 0.7e0],
            [0.0e0, -43.6e0, 2.1e0, 0.0e0],
            [0.0e0, -25.4e0, 1.2e0, 0.0e0],
            [2.0e0, 40.9e0, -2.0e0, 0.0e0],
            [-2.1e0, 26.1e0, 0.6e0, 0.0e0],
            [22.6e0, -3.2e0, -0.5e0, -0.5e0],
            [-7.6e0, 24.9e0, -0.4e0, -0.2e0],
            [-6.2e0, 34.9e0, 1.7e0, 0.3e0],
            [2.0e0, 17.4e0, -0.4e0, 0.1e0],
            [-3.9e0, 20.5e0, 2.4e0, 0.6e0],
        ])

        # Nutation series: obliquity
        EPS = np.array([
            [9205365.8e0, -1506.2e0, 885.7e0, -0.2e0],
            [573095.9e0, -570.2e0, -305.0e0, -0.3e0],
            [97845.5e0, 147.8e0, -48.8e0, -0.2e0],
            [-89753.6e0, 28.0e0, 46.9e0, 0.0e0],
            [7406.7e0, -327.1e0, -18.2e0, 0.8e0],
            [22442.3e0, -22.3e0, -67.6e0, 0.0e0],
            [-683.6e0, 46.8e0, 0.0e0, 0.0e0],
            [20070.7e0, 36.0e0, 1.6e0, 0.0e0],
            [12893.8e0, 39.5e0, -6.2e0, 0.0e0],
            [-9593.2e0, 14.4e0, 30.2e0, -0.1e0],
            [-6899.5e0, 4.8e0, -0.6e0, 0.0e0],
            [-5332.5e0, -0.1e0, 2.7e0, 0.0e0],
            [-125.2e0, 10.5e0, 0.0e0, 0.0e0],
            [-3323.4e0, -0.9e0, -0.3e0, 0.0e0],
            [3142.3e0, 8.9e0, 0.3e0, 0.0e0],
            [2552.5e0, 7.3e0, -1.2e0, 0.0e0],
            [2634.4e0, 8.8e0, 0.2e0, 0.0e0],
            [-2424.4e0, 1.6e0, -0.4e0, 0.0e0],
            [-123.3e0, 3.9e0, 0.0e0, 0.0e0],
            [1642.4e0, 7.3e0, -0.8e0, 0.0e0],
            [47.9e0, 3.2e0, 0.0e0, 0.0e0],
            [1321.2e0, 6.2e0, -0.6e0, 0.0e0],
            [-1234.1e0, -0.3e0, 0.6e0, 0.0e0],
            [-1076.5e0, -0.3e0, 0.0e0, 0.0e0],
            [-61.6e0, 1.8e0, 0.0e0, 0.0e0],
            [-55.4e0, 1.6e0, 0.0e0, 0.0e0],
            [856.9e0, -4.9e0, -2.1e0, 0.0e0],
            [-800.7e0, -0.1e0, 0.0e0, 0.0e0],
            [685.1e0, -0.6e0, -3.8e0, 0.0e0],
            [-16.9e0, -1.5e0, 0.0e0, 0.0e0],
            [695.7e0, 1.8e0, 0.0e0, 0.0e0],
            [642.2e0, -2.6e0, -1.6e0, 0.0e0],
            [13.3e0, 1.1e0, -0.1e0, 0.0e0],
            [521.9e0, 1.6e0, 0.0e0, 0.0e0],
            [325.8e0, 2.0e0, -0.1e0, 0.0e0],
            [-325.1e0, -0.5e0, 0.9e0, 0.0e0],
            [10.1e0, 0.3e0, 0.0e0, 0.0e0],
            [334.5e0, 1.6e0, 0.0e0, 0.0e0],
            [307.1e0, 0.4e0, -0.9e0, 0.0e0],
            [327.2e0, 0.5e0, 0.0e0, 0.0e0],
            [-304.6e0, -0.1e0, 0.0e0, 0.0e0],
            [304.0e0, 0.6e0, 0.0e0, 0.0e0],
            [-276.8e0, -0.5e0, 0.1e0, 0.0e0],
            [268.9e0, 1.3e0, 0.0e0, 0.0e0],
            [271.8e0, 1.1e0, 0.0e0, 0.0e0],
            [271.5e0, -0.4e0, -0.8e0, 0.0e0],
            [-5.2e0, 0.5e0, 0.0e0, 0.0e0],
            [-220.5e0, 0.1e0, 0.0e0, 0.0e0],
            [-20.1e0, 0.3e0, 0.0e0, 0.0e0],
            [-191.0e0, 0.1e0, 0.5e0, 0.0e0],
            [-4.1e0, 0.3e0, 0.0e0, 0.0e0],
            [130.6e0, -0.1e0, 0.0e0, 0.0e0],
            [3.0e0, 0.3e0, 0.0e0, 0.0e0],
            [122.9e0, 0.8e0, 0.0e0, 0.0e0],
            [3.7e0, -0.3e0, 0.0e0, 0.0e0],
            [123.1e0, 0.4e0, -0.3e0, 0.0e0],
            [-52.7e0, 15.3e0, 0.0e0, 0.0e0],
            [120.7e0, 0.3e0, -0.3e0, 0.0e0],
            [4.0e0, -0.3e0, 0.0e0, 0.0e0],
            [126.5e0, 0.5e0, 0.0e0, 0.0e0],
            [112.7e0, 0.5e0, -0.3e0, 0.0e0],
            [-106.1e0, -0.3e0, 0.3e0, 0.0e0],
            [-112.9e0, -0.2e0, 0.0e0, 0.0e0],
            [3.6e0, -0.2e0, 0.0e0, 0.0e0],
            [107.4e0, 0.3e0, 0.0e0, 0.0e0],
            [-10.9e0, 0.2e0, 0.0e0, 0.0e0],
            [-0.9e0, 0.0e0, 0.0e0, 0.0e0],
            [85.4e0, 0.0e0, 0.0e0, 0.0e0],
            [0.0e0, -88.8e0, 0.0e0, 0.0e0],
            [-71.0e0, -0.2e0, 0.0e0, 0.0e0],
            [-70.3e0, 0.0e0, 0.0e0, 0.0e0],
            [64.5e0, 0.4e0, 0.0e0, 0.0e0],
            [69.8e0, 0.0e0, 0.0e0, 0.0e0],
            [66.1e0, 0.4e0, 0.0e0, 0.0e0],
            [-61.0e0, -0.2e0, 0.0e0, 0.0e0],
            [-59.5e0, -0.1e0, 0.0e0, 0.0e0],
            [-55.6e0, 0.0e0, 0.2e0, 0.0e0],
            [51.7e0, 0.2e0, 0.0e0, 0.0e0],
            [-49.0e0, -0.1e0, 0.0e0, 0.0e0],
            [-52.7e0, -0.1e0, 0.0e0, 0.0e0],
            [-49.6e0, 1.4e0, 0.0e0, 0.0e0],
            [46.3e0, 0.4e0, 0.0e0, 0.0e0],
            [49.6e0, 0.1e0, 0.0e0, 0.0e0],
            [-5.1e0, 0.1e0, 0.0e0, 0.0e0],
            [-44.0e0, -0.1e0, 0.0e0, 0.0e0],
            [-39.9e0, -0.1e0, 0.0e0, 0.0e0],
            [-39.5e0, -0.1e0, 0.0e0, 0.0e0],
            [-3.9e0, 0.1e0, 0.0e0, 0.0e0],
            [-42.1e0, -0.1e0, 0.0e0, 0.0e0],
            [-17.2e0, 0.1e0, 0.0e0, 0.0e0],
            [-2.3e0, 0.1e0, 0.0e0, 0.0e0],
            [-39.2e0, 0.0e0, 0.0e0, 0.0e0],
            [-38.4e0, 0.1e0, 0.0e0, 0.0e0],
            [36.8e0, 0.2e0, 0.0e0, 0.0e0],
            [34.6e0, 0.1e0, 0.0e0, 0.0e0],
            [-32.7e0, 0.3e0, 0.0e0, 0.0e0],
            [30.4e0, 0.0e0, 0.0e0, 0.0e0],
            [0.4e0, 0.1e0, 0.0e0, 0.0e0],
            [29.3e0, 0.2e0, 0.0e0, 0.0e0],
            [31.6e0, 0.1e0, 0.0e0, 0.0e0],
            [0.8e0, -0.1e0, 0.0e0, 0.0e0],
            [-27.9e0, 0.0e0, 0.0e0, 0.0e0],
            [2.9e0, 0.0e0, 0.0e0, 0.0e0],
            [-25.3e0, 0.0e0, 0.0e0, 0.0e0],
            [25.0e0, 0.1e0, 0.0e0, 0.0e0],
            [27.5e0, 0.1e0, 0.0e0, 0.0e0],
            [-24.4e0, -0.1e0, 0.0e0, 0.0e0],
            [24.9e0, 0.2e0, 0.0e0, 0.0e0],
            [-22.8e0, -0.1e0, 0.0e0, 0.0e0],
            [0.9e0, -0.1e0, 0.0e0, 0.0e0],
            [24.4e0, 0.1e0, 0.0e0, 0.0e0],
            [23.9e0, 0.1e0, 0.0e0, 0.0e0],
            [22.5e0, 0.1e0, 0.0e0, 0.0e0],
            [20.8e0, 0.1e0, 0.0e0, 0.0e0],
            [20.1e0, 0.0e0, 0.0e0, 0.0e0],
            [21.5e0, 0.1e0, 0.0e0, 0.0e0],
            [-20.0e0, 0.0e0, 0.0e0, 0.0e0],
            [1.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-0.2e0, -0.1e0, 0.0e0, 0.0e0],
            [19.0e0, 0.0e0, -0.1e0, 0.0e0],
            [20.5e0, 0.0e0, 0.0e0, 0.0e0],
            [-2.0e0, 0.0e0, 0.0e0, 0.0e0],
            [-17.6e0, -0.1e0, 0.0e0, 0.0e0],
            [19.0e0, 0.0e0, 0.0e0, 0.0e0],
            [-2.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-18.4e0, -0.1e0, 0.0e0, 0.0e0],
            [17.1e0, 0.0e0, 0.0e0, 0.0e0],
            [0.4e0, 0.0e0, 0.0e0, 0.0e0],
            [18.4e0, 0.1e0, 0.0e0, 0.0e0],
            [0.0e0, 17.4e0, 0.0e0, 0.0e0],
            [-0.6e0, 0.0e0, 0.0e0, 0.0e0],
            [-15.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-16.8e0, -0.1e0, 0.0e0, 0.0e0],
            [16.3e0, 0.0e0, 0.0e0, 0.0e0],
            [-2.0e0, 0.0e0, 0.0e0, 0.0e0],
            [-1.5e0, 0.0e0, 0.0e0, 0.0e0],
            [-14.3e0, -0.1e0, 0.0e0, 0.0e0],
            [14.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-13.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-14.3e0, -0.1e0, 0.0e0, 0.0e0],
            [-13.7e0, 0.0e0, 0.0e0, 0.0e0],
            [13.1e0, 0.1e0, 0.0e0, 0.0e0],
            [-1.7e0, 0.0e0, 0.0e0, 0.0e0],
            [-12.8e0, 0.0e0, 0.0e0, 0.0e0],
            [0.0e0, -14.4e0, 0.0e0, 0.0e0],
            [12.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-12.0e0, 0.0e0, 0.0e0, 0.0e0],
            [-0.8e0, 0.0e0, 0.0e0, 0.0e0],
            [10.9e0, 0.1e0, 0.0e0, 0.0e0],
            [-10.8e0, 0.0e0, 0.0e0, 0.0e0],
            [10.5e0, 0.0e0, 0.0e0, 0.0e0],
            [-10.4e0, 0.0e0, 0.0e0, 0.0e0],
            [-11.2e0, 0.0e0, 0.0e0, 0.0e0],
            [10.5e0, 0.1e0, 0.0e0, 0.0e0],
            [-1.4e0, 0.0e0, 0.0e0, 0.0e0],
            [0.0e0, 0.1e0, 0.0e0, 0.0e0],
            [0.7e0, 0.0e0, 0.0e0, 0.0e0],
            [-10.3e0, 0.0e0, 0.0e0, 0.0e0],
            [-10.0e0, 0.0e0, 0.0e0, 0.0e0],
            [9.6e0, 0.0e0, 0.0e0, 0.0e0],
            [9.4e0, 0.1e0, 0.0e0, 0.0e0],
            [0.6e0, 0.0e0, 0.0e0, 0.0e0],
            [-87.7e0, 4.4e0, -0.4e0, -6.3e0],
            [46.3e0, 22.4e0, 0.5e0, -2.4e0],
            [15.6e0, -3.4e0, 0.1e0, 0.4e0],
            [5.2e0, 5.8e0, 0.2e0, -0.1e0],
            [-30.1e0, 26.9e0, 0.7e0, 0.0e0],
            [23.2e0, -0.5e0, 0.0e0, 0.6e0],
            [1.0e0, 23.2e0, 3.4e0, 0.0e0],
            [-12.2e0, -4.3e0, 0.0e0, 0.0e0],
            [-2.1e0, -3.7e0, -0.2e0, 0.1e0],
            [-18.6e0, -3.8e0, -0.4e0, 1.8e0],
            [5.5e0, -18.7e0, -1.8e0, -0.5e0],
            [-5.5e0, -18.7e0, 1.8e0, -0.5e0],
            [18.4e0, -3.6e0, 0.3e0, 0.9e0],
            [-0.6e0, 1.3e0, 0.0e0, 0.0e0],
            [-5.6e0, -19.5e0, 1.9e0, 0.0e0],
            [5.5e0, -19.1e0, -1.9e0, 0.0e0],
            [-17.3e0, -0.8e0, 0.0e0, 0.9e0],
            [-3.2e0, -8.3e0, -0.8e0, 0.3e0],
            [-0.1e0, 0.0e0, 0.0e0, 0.0e0],
            [-5.4e0, 7.8e0, -0.3e0, 0.0e0],
            [-14.8e0, 1.4e0, 0.0e0, 0.3e0],
            [-3.8e0, 0.4e0, 0.0e0, -0.2e0],
            [12.6e0, 3.2e0, 0.5e0, -1.5e0],
            [0.1e0, 0.0e0, 0.0e0, 0.0e0],
            [-13.6e0, 2.4e0, -0.1e0, 0.0e0],
            [0.9e0, 1.2e0, 0.0e0, 0.0e0],
            [-11.9e0, -0.5e0, 0.0e0, 0.3e0],
            [0.4e0, 12.0e0, 0.3e0, -0.2e0],
            [8.3e0, 6.1e0, -0.1e0, 0.1e0],
            [0.0e0, 0.0e0, 0.0e0, 0.0e0],
            [0.4e0, -10.8e0, 0.3e0, 0.0e0],
            [9.6e0, 2.2e0, 0.3e0, -1.2e0],

        ])

        # Interval between fundamental epoch J2000.0 and given epoch (JC).
        T = (date - cls.DJM0) / cls.DJC

        # Mean anomaly of the Moon.
        EL = 134.96340251e0 * cls.DD2R + np.mod(T * (1717915923.2178e0 +
                                                     T * (31.8792e0 +
                                                          T * (0.051635e0 +
                                                               T * (- 0.00024470e0)))), cls.TURNAS) * cls.AS2R

        # Mean anomaly of the Sun.
        ELP = 357.52910918e0 * cls.DD2R + np.mod(T * (129596581.0481e0 +
                                                      T * (- 0.5532e0 +
                                                           T * (0.000136e0 +
                                                                T * (- 0.00001149e0)))), cls.TURNAS) * cls.AS2R

        # Mean argument of the latitude of the Moon.
        F = 93.27209062e0 * cls.DD2R + np.mod(T * (1739527262.8478e0 +
                                                   T * (- 12.7512e0 +
                                                        T * (- 0.001037e0 +
                                                             T * 0.00000417e0))), cls.TURNAS) * cls.AS2R

        # Mean elongation of the Moon from the Sun.
        D = 297.85019547e0 * cls.DD2R + np.mod(T * (1602961601.2090e0 +
                                                    T * (- 6.3706e0 +
                                                         T * (0.006539e0 +
                                                              T * (- 0.00003169e0)))), cls.TURNAS) * cls.AS2R

        # Mean longitude of the ascending node of the Moon.
        OM = 125.04455501e0 * cls.DD2R + np.mod(T * (- 6962890.5431e0 +
                                                     T * (7.4722e0 +
                                                          T * (0.007702e0 +
                                                               T * (- 0.00005939e0)))), cls.TURNAS) * cls.AS2R

        # Mean longitude of Venus.
        VE = 181.97980085e0 * cls.DD2R + np.mod(210664136.433548e0 * T, cls.TURNAS) * cls.AS2R

        # Mean longitude of Mars.
        MA = 355.43299958e0 * cls.DD2R + np.mod(68905077.493988e0 * T, cls.TURNAS) * cls.AS2R

        # Mean longitude of Jupiter.
        JU = 34.35151874e0 * cls.DD2R + np.mod(10925660.377991e0 * T, cls.TURNAS) * cls.AS2R

        # Mean longitude of Saturn.
        SA = 50.07744430e0 * cls.DD2R + np.mod(4399609.855732e0 * T, cls.TURNAS) * cls.AS2R

        # Geodesic nutation (Fukushima 1991) in microarcsec.
        DP = -153.1e0 * np.sin(ELP) - 1.9e0 * np.sin(2e0 * ELP)
        DE = 0e0

        # Shirai & Fukushima (2001) nutation series.

        na_flipped = np.flip(NA, 0)
        psi_flipped = np.flip(PSI, 0)
        eps_flipped = np.flip(EPS, 0)
        coeff_list = np.array([EL, ELP, F, D, OM, VE, MA, JU, SA])
        for na_row, psi_row, eps_row in zip(na_flipped, psi_flipped, eps_flipped):
            THETA = np.dot(na_row, coeff_list)

            C = np.cos(THETA)
            S = np.sin(THETA)
            DP = DP + (psi_row[0] + psi_row[2] * T) * C + (psi_row[1] + psi_row[3] * T) * S
            DE = DE + (eps_row[0] + eps_row[2] * T) * C + (eps_row[1] + eps_row[3] * T) * S

        # Change of units, and addition of the precession correction.
        dpsi = (DP * 1e-6 - 0.042888e0 - 0.29856e0 * T) * cls.AS2R
        deps = (DE * 1e-6 - 0.005171e0 - 0.02408e0 * T) * cls.AS2R

        # Mean obliquity of date (Simon et al. 1994).
        eps0 = (84381.412e0 + (-46.80927e0 +
                               (-0.000152e0 +
                                (0.0019989e0 +
                                 (-0.00000051e0 +
                                  (-0.000000025e0) * T) * T) * T) * T) * T) * cls.AS2R

        return dpsi, deps, eps0

    @classmethod
    def vdv(cls, in_vec_a, in_vec_b):
        """+
        - - - -
        V D V
        - - - -

        Scalar product of two 3-vectors  (single precision)

        Given:
        VA      real(3)     first vector
        VB      real(3)     second vector

        The result is the scalar product VA.VB (single precision)
        -"""
    
        return np.dot(in_vec_a, in_vec_b)
    
    @classmethod
    def ds2tp(cls, ra, dec, raz, decz):
        """
        - - - - - -
        D S 2 T P
        - - - - - -

        Projection of spherical coordinates onto tangent plane:
        "gnomonic" projection - "standard coordinates" (double precision)

        Given:
        RA,DEC      dp   spherical coordinates of point to be projected
        RAZ,DECZ    dp   spherical coordinates of tangent point

        Returned:
        XI,ETA      dp   rectangular coordinates on tangent plane
        J           int  status:   0 = OK, star on tangent plane
        1 = error, star too far from axis
        2 = error, antistar on tangent plane
        3 = error, antistar too far from axis
        Trig functions
        """
        Sdecz = np.sin(decz)
        Sdec = np.sin(dec)
        Cdecz = np.cos(decz)
        Cdec = np.cos(dec)
        RADIF = ra - raz
        SRADIF = np.sin(RADIF)
        CRADIF = np.cos(RADIF)

        # Reciprocal of star vector length to tangent plane
        DENOM = Sdec*Sdecz+Cdec*Cdecz+CRADIF

        # Handle vectors too far from axis
        if DENOM > 1e-6:
            j = 0
        elif DENOM >= 0e0:
            j = 1
            DENOM = 1e-6
        elif DENOM > -1e-6:
            j = 2
            DENOM = -1e-6
        else:
            j = 3

        # Compute tangent plane coordinates (even in dubious cases)
        xi = Cdec*SRADIF/DENOM
        eta = (Sdec*Cdecz-Cdec*Sdecz-CRADIF)/DENOM
        return xi, eta, j
