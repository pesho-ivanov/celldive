import os
import shutil
import sys 
import argparse

DEFAULT_CONF = os.path.join(os.path.expanduser("~"), '.bioinf.conf')

#import site;
#print (site.getsitepackages())

# TODO: http://stackoverflow.com/questions/31949760/how-to-limit-python-traceback-to-specific-files
# Don't print traceback on exceptions
#sys.tracebacklimit = 1

def get_config(conf_file):
    def toDict(category):
        return {key: value for (key, value) in category.items()}

    config = configparser.ConfigParser()
    config.read(conf_file)
    libs = config['libs']
    data = config['data']
    print(libs.keys())
    return (toDict(libs), toDict(data))

def run_cmd(cmd):
    print("Envoking: ", '\n'.join(cmd))
    subprocess.call(cmd)

#def to_full_path(path):
#    return path
#    #return os.path.abspath(os.path.expanduser(path))

# hack from https://gist.github.com/brantfaircloth/1443543
def VerifyFile(extensions={}, should_exist=True):
    """Checks """
    class Act(argparse.Action):
        def __call__(self, parser, namespace, fname, option_string=None):
            if should_exist and not os.path.exists(fname):
                parser.error("{0} is not a file".format(fname))
            ext = os.path.splitext(fname)[1][1:]
            if ext.lower() not in map(ext.lower(), extensions):
                option_string = '({})'.format(option_string) if option_string else ''
                parser.error("file doesn't end with one of {}{}".format(extensions, option_string))
            else:
                setattr(namespace, self.dest, fname)
    return Act

# from http://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input
def query_yes_no(question, default=None):
    """Ask a yes/no question via input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes", "no" or None (meaning an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def create_dir(dirname):
    """Create a directory or clean it with the permission of the user."""
    dirname = to_full_path(dirname)
    if os.path.isdir(dirname):
        files = os.listdir(dirname)
        if len(files) > 0:
            question = "{0} already exists: Do you want to clear it?".format(dirname)
            if not query_yes_no(question, "no"): \
                raise ValueError("Clear the output dir yourself or choose another one") 
               
            print('\n'.join([ ' > ' + d for d in files ]))
            if not query_yes_no("Are you sure you want to delete all this?", "no"):
                raise IOError("Clear the output dir yourself or choose another one")
            print(dirname + " successfully cleaned.")
        shutil.rmtree(dirname)

    os.makedirs(dirname)
    if dirname[-1] != '/':
        dirname += '/'
    return dirname

def verify_file_exists(filename):
    """Checks if a path is an actual file"""
    if not os.path.exists(filename):
        raise IOError("{0} is not a file".format(filename))
    return filename 

def match_file_type(filename, ext, gz=False):
    if gz:
        for e in list(ext):
            ext.append(e + '.gz')
    if not filename.endswith(tuple(ext)):
        raise ValueError("The file " + filename + " doens't have any of the types " + ','.join(ext))
    return verify_file_exists(filename)

def match_fasta(filename, gz=False):
    return match_file_type(filename, ['.fasta', '.fa'], gz)

def match_fastq(filename, gz=True):
    return match_file_type(filename, ['.fastq', '.fq'], gz)

def match_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        raise NotADirectoryError("{0} is not a directory".format(dirname))
    return dirname

def is_gz(fn):
    assert type(fn) is str
    return fn.lower().endswith('.gz')

def are_gz(files):
    assert type(fn) is list
    is_gz = [fn.lower().endswith('.gz') for fn in files]
    if len(set(is_gz)) > 1:
       raise IOError("Only some files are .gz")
    print ('GZ: ', is_gz[0])
    return is_gz[0]
