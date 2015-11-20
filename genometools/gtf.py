"""Functions for GTF files."""

def parse_attributes(s):
    """ Parses the ``attribute`` string of a GFF/GTF annotation.

    Parameters
    ----------
    s : str
        The attribute string.

    Returns
    -------
    Dict
        A dictionary containing attribute name/value pairs.

    Notes
    -----
    The ``attribute`` string is the 9th field of each annotation (row),
    as described in the
    `GTF format specification <http://mblab.wustl.edu/GTF22.html>`_.

    """

    # use regular expression with negative lookbehind to make sure we don't
    # split on escaped semicolons ("\;")
    attr_sep = re.compile(r"(?<!\\)\s*;\s*")
    attr = {}
    atts = attr_sep.split(s)
    for a in atts:
        #print a
        kv = a.split(' ')
        if len(kv) == 2:
            k,v = kv
            v = v.strip('"')
            attr[k] = v
    return attr 


