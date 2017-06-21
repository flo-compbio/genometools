import ftputil

sra_protocol = 'ftp://'
sra_host = 'ftp-trace.ncbi.nlm.nih.gov'
sra_user = 'anonymous'
sra_password = 'anonymous'

def get_experiment_urls(exp):
    exp_urls = []
    with ftputil.FTPHost(sra_host, sra_user, sra_password) as ftp_host:
        download_paths = []
        exp_dir = '/sra/sra-instant/reads/ByExp/sra/SRX/%s/%s/' \
                %(exp[:6],exp)
        ftp_host.chdir(exp_dir)
        run_folders = ftp_host.listdir(ftp_host.curdir)

        # compile a list of all files
        for folder in run_folders:
            files = ftp_host.listdir(folder)
            assert len(files) == 1
            for f in files:
                path = exp_dir + folder + '/' + f
                exp_urls.append(path)
    return exp_urls


def get_project_urls(project):
    """Get the URLs for all runs from a given project.
    
    TODO: docstring"""
    urls = []
    with ftputil.FTPHost(sra_host, sra_user, sra_password) as ftp_host:
        download_paths = []
        exp_dir = '/sra/sra-instant/reads/ByStudy/sra/SRP/%s/%s/' \
                %(project[:6], project)
        ftp_host.chdir(exp_dir)
        run_folders = ftp_host.listdir(ftp_host.curdir)

        # compile a list of all files
        for folder in run_folders:
            files = ftp_host.listdir(folder)
            assert len(files) == 1
            for f in files:
                path = exp_dir + folder + '/' + f
                urls.append(path)
    return urls


def get_run_url(run):
    """Returns the FTP URL for a given SRA run accession ID."""
    return ('ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra'
            '/SRR/%s/%s/%s.sra' % (run[:6], run, run))


def get_file_sizes(files):
    sizes = []
    with ftputil.FTPHost(sra_host, sra_user, sra_password) as ftp_host:
        for f in files:
            stat = ftp_host.stat(f)
            sizes.append(stat.st_size)
    return sizes
