"""Miscellaneous functions used by several classes"""
def isoformat_digits(time, digits):
    """
    return datetime as isoformat with specified digits after decimal
    
    :param digits: 0-6
    """
    if not time:
        return 'None'
    s = time.strftime('%Y-%m-%dT%H:%M:%S')
    digits = int(digits)
    if digits <= 0:
        return s
    if digits > 6:
        digits = 6
    fmt='.{:0' + str(digits) + 'd}'
    s += fmt.format(int(time.microsecond * 10**(digits-6)))
    return s
