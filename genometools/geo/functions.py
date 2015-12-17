def get_series_matrix_url(a):
    return 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz' \
            %(a[:5] + 'nnn', a, a)

def get_raw_data_url(a):
    return 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s_RAW.tar' \
            %(a[:5] + 'nnn', a, a)
