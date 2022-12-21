"""
A date time parse mostly written for parsing the PATRIC metadata

Note: for a functional version, please see the jupyter notebook "RAST sources.ipynb" in the PhispyAnalysis github repo.
"""
import sys

from dateutil.parser import parse, ParserError
import pytz

import re

# Convert the collection date to a number. We take todays date, and then convert relative to that
# it allows dates pre-epoch

yr = re.compile("^[~]*(\\d{4})['s]*$")
yrrange = re.compile("^(\\d{4})\\s*.\\s*(\\d{4})$")
myd = re.compile('^\\d{1,4}/\\d{1,2}/\\d{1,4}$')
damoyrrange = re.compile('^(\\d+\\s+\\w+\\s+\\d{4})\\s*.{1,3}?\\s*(\\d+\\s+\\w+\\s+\\d{4})$')
modacyrrange = re.compile('^(\\w+\\s+\\d+,*\\s+\\d{4})\\s*.{1,3}?\\s*(\\w+\\s+\\d+,*\\s+\\d{4})$')
moyrrange = re.compile('^(\\w+\\s+\\d{4})\\s*.{1,3}?\\s*(\\w+\\s+\\d{4})$')
year42 = re.compile('^(\\d{4})-\\d{2}$')
lem = re.compile("^late\\s*|^early\\s*|^mid\\s*|^prior to\\s*|^before\\s*|^pre-", re.IGNORECASE)
splitseen = set()


def _try_parsing(test_date, adate):
    """
    Attempt to parse a date, and catch an error.

    If we fail, we return None, otherwise we return the years since now()
    """
    try:
        dt = parse(test_date)
    except:
        return None

    dt = dt.replace(tzinfo=pytz.UTC)

    if dt < adate:
        tdelt = adate - dt
        seconds = -1 * ((tdelt.days * 86400) + tdelt.seconds)
    else:
        tdelt = dt - adate
        seconds = (tdelt.days * 86400) + tdelt.seconds
    # convert seconds to years
    return seconds / 31557600


def date_to_seconds(test_date, epoch="24/02/1977"):
    """
    Convert the date to years and fractions.

    We try several times, and clean it up as we go along.
    :param test_date: the date to try
    :param epoch: the date to set as year 0 (things before will be negative, things after positive)
    """
    if not test_date:
        # may need to check for np.nan here
        return None

    try:
        adate = parse(epoch)
    except ParserError as e:
        sys.stderr.write(f"FATAL: {epoch} is not a valid date ({e})\n")
        sys.exit(-1)
    adate = adate.replace(tzinfo=pytz.UTC)

    # can we parse this date? If so, lets do it and return the value
    attempt = _try_parsing(test_date, adate)
    if attempt:
        return attempt
    orix = test_date

    if test_date.lower() in ['restricted access', 'none', 'not collected', 'not applicable',
                             'not available: not collected', 'unspecified']:
        return None

    # a few one off cases that are just easier to fix
    if 'May 2015-Nov 2015' == test_date:
        test_date = 'May 2015'

    if '1954-65' == test_date:
        test_date = '01 January 1954'

    if 'Jul-00' == test_date:
        test_date = 'Jul-2000'

    if '2015_9' == test_date:
        test_date = 'Sep-2015'

    if '31-Mac-2013' == test_date:
        test_date = '31-May-2013'

    if '2010-0916' == test_date:
        test_date = '16 Sep 2010'

    test_date = lem.sub('', test_date)

    if '_' in test_date:
        test_date = test_date.replace('_', '-')

    test_date = test_date.replace(' or earlier', '')
    test_date = test_date.replace('collected in the ', '')

    # some regular expressions of variants of day month year - day month year ranges. We choose 1
    m = yrrange.match(test_date)
    if m:
        test_date = '01 January ' + m.groups()[1]

    m = yr.match(test_date)
    if m:
        test_date = '01 January ' + m.groups()[0]

    m = year42.match(test_date)
    if m:
        test_date = '01 January ' + m.groups()[0]

    m = modacyrrange.match(test_date)
    if m:
        test_date = m.groups()[1]

    m = damoyrrange.match(test_date)
    if m:
        test_date = m.groups()[1]

    m = moyrrange.match(test_date)
    if m:
        test_date = '01 ' + m.groups()[1]

    # can we parse this date? If so, lets do it and return the value
    attempt = _try_parsing(test_date, adate)
    if attempt:
        return attempt

    if '/' in test_date:
        if test_date not in splitseen:
            # sys.stderr.write(f"Splitting {x}\n")
            splitseen.add(test_date)
        p = test_date.split('/')
        test_date = p[1]

    # can we parse this date? If so, lets do it and return the value
    attempt = _try_parsing(test_date, adate)
    if attempt:
        return attempt

    if test_date.endswith('-00'):
        test_date = test_date.replace('-00', '-2000')

    # can we parse this date? If so, lets do it and return the value
    attempt = _try_parsing(test_date, adate)
    if attempt:
        return attempt

    sys.stderr.write(f"can't parse |{test_date}| from |{orix}|\n")

    return None
