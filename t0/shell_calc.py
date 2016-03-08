from __future__ import division
from collections import deque
from os import getcwd, path, walk, mkdir, link, chdir
from subprocess import call
import re
import glob

# CONSTANTS
# .ans file
SEP = '--------------------------------------------------'
LINES = ['%s',
         '%s,   %d',
         '%s',
         '%s',
         '%s',
         ' %d',
         ' %d',
         ' %.1f, %.1f, %.1f',
         ' %d',
         '%s',
         '%s']
NUM_PROTONS = 8

# directories
D = getcwd()
SOURCES = path.join(D, 'sources')
RESULTS = path.join(D, 'results')
DIRNAME_USDB = 'usdb'
MODEL_SPACE = 'imsrg'
MODEL_SPACE_USDB = 'sd'

# file parsing
MASS_REGEX = 'A\d+'
FILENAME_SPLIT = '_'
EXT_REGEX = '.*int'
ANS_REGEX = '.*ans'
BAT_REGEX = MASS_REGEX + '\.bat'

MASS_RANGE = range(17, 29)


def do_all_calculations(arange=MASS_RANGE, zrange=list([NUM_PROTONS]),
                        **kwargs):
    for z in zrange:
        make_results_dir(num_protons=z, **kwargs)
        make_usdb_dir(mass_range=arange, num_protons=z, **kwargs)
    return do_calculations(**kwargs)


def do_calculations(d=D, results_dir=RESULTS,
                    ans_regex=ANS_REGEX, bat_regex=BAT_REGEX,
                    force=False):
    for root, dirs, files in walk(results_dir):
        ans = _get_file(files, ans_regex)
        bat = _get_file(files, bat_regex)
        if ans is not None:  # There is a *.ans file
            if bat is None or force:
                chdir(root)
                call(['shell', '%s' % ans])
    chdir(d)
    for root, dirs, files in walk(results_dir):
        bat = _get_file(files, bat_regex)
        if bat is not None:
            if not _calc_has_been_done(root) or force is True:
                chdir(root)
                try:
                    call(['source', '%s' % bat])
                except OSError:
                    call(['source %s' % bat], shell=True)
    chdir(d)
    return 1


def _calc_has_been_done(dirpath):
    return len(files_with_ext_in_directory(dirpath, '.lpt')) > 1


def files_with_ext_in_directory(directory, extension):
    """Returns a list of the filenames of all the files_INT in the given
    directory with the given extension
    :param extension:
    :param directory: """
    return list(glob.glob(path.join(directory, '*' + extension)))


def _get_file(list_of_file, regex=ANS_REGEX):
    for f in list_of_file:
        if re.match(regex, f) is not None:
            return f
    else:
        return None


def make_usdb_dir(mass_range=MASS_RANGE, d=D, results_dir=RESULTS,
                  model_space=MODEL_SPACE_USDB,
                  num_protons=NUM_PROTONS,
                  usdb_dirname=DIRNAME_USDB,
                  force=False):
    results_subdir = path.join(results_dir, 'Z%d' % num_protons)
    usdb_dirpath = path.join(results_subdir, usdb_dirname)
    if not path.exists(usdb_dirpath):
        mkdir(usdb_dirpath)
    for mass_num in mass_range:
        dirname = path.join(usdb_dirpath, 'A%d' % mass_num)
        if not path.exists(dirname):
            mkdir(dirname)
        # link .sp file
        sp_filename = '%s.sp' % model_space
        sp_file_path = path.join(dirname, sp_filename)
        if not path.exists(sp_file_path):
            link(path.join(d, sp_filename), sp_file_path)
        # create .ans file
        ans_filename = 'A%d.ans' % mass_num
        ans_file_path = path.join(dirname, ans_filename)
        if not path.exists(ans_file_path) or force is True:
            if mass_num % 2 == 0:  # even
                make_ans_file(file_path=ans_file_path,
                              sp_file=model_space,
                              num_nucleons=mass_num,
                              interaction_name='usdb',
                              num_protons=num_protons)
            else:
                make_ans_file(file_path=ans_file_path,
                              sp_file=model_space,
                              num_nucleons=mass_num,
                              interaction_name='usdb',
                              num_protons=num_protons,
                              j_min=0.5, j_max=3.5, j_del=1.0)


def make_results_dir(d=D, sources_dir=SOURCES, results_dir=RESULTS,
                     model_space=MODEL_SPACE,
                     file_ext_regex=EXT_REGEX,
                     num_protons=NUM_PROTONS,
                     force=False):
    """Copy all of the directories from the sources_dir into the results_dir
    recursively, but when encountering a *.int file, make a directory for the
    file (according to its name) and copy the file into the directory with a
    short name. Also, for each directory to which a *.int file is copied the
    given model space is linked and a *.ans file is generated.
    :param force:
    :param num_protons:
    :param file_ext_regex:
    :param model_space:
    :param results_dir:
    :param sources_dir:
    :param d:
    """
    results_subdir = path.join(results_dir, 'Z%d' % num_protons)
    if not path.exists(sources_dir):
        raise SourcesDirDoesNotExistException()
    if not path.exists(results_subdir):
        mkdir(results_subdir)

    todo_sources = deque([sources_dir])
    todo_results = deque([results_subdir])

    while len(todo_sources) > 0:
        cwd_sources = todo_sources.popleft()
        cwd_results = todo_results.popleft()
        root, dirs, files = walk(cwd_sources).next()
        for dd in dirs:
            next_sources = path.join(cwd_sources, dd)
            next_results = path.join(cwd_results, dd)
            if not path.exists(next_results):
                mkdir(next_results)
            todo_sources.append(next_sources)
            todo_results.append(next_results)
        for ff in filter(lambda f: re.match(file_ext_regex, f) is not None,
                         files):
            new_dir = path.join(cwd_results, _fname_without_extension(ff))
            if not path.exists(new_dir):
                mkdir(new_dir)
            mass_num = mass_number_from_filename(ff)
            # link .int file
            interaction_name = 'A%d' % mass_num
            interaction_file_path = path.join(new_dir,
                                              interaction_name + '.int')
            if not path.exists(interaction_file_path):
                link(path.join(cwd_sources, ff), interaction_file_path)
            # link .sd file
            sp_filename = '%s.sp' % model_space
            sp_file_path = path.join(new_dir, sp_filename)
            if not path.exists(sp_file_path):
                link(path.join(d, sp_filename), sp_file_path)
            # create .ans file
            ans_filename = 'A%d.ans' % mass_num
            ans_file_path = path.join(new_dir, ans_filename)
            if not path.exists(ans_file_path) or force is True:
                if mass_num % 2 == 0:  # even
                    make_ans_file(file_path=ans_file_path,
                                  option='lpe', neig=0,
                                  sp_file=model_space,
                                  restriction='n',
                                  interaction_name=interaction_name,
                                  num_protons=num_protons,
                                  num_nucleons=mass_num,
                                  j_min=0.0, j_max=4.0, j_del=1.0,
                                  parity=0)
                else:
                    make_ans_file(file_path=ans_file_path,
                                  option='lpe', neig=0,
                                  sp_file=model_space,
                                  restriction='n',
                                  interaction_name=interaction_name,
                                  num_protons=num_protons,
                                  num_nucleons=mass_num,
                                  j_min=0.5, j_max=3.5, j_del=1.0,
                                  parity=0)


def mass_number_from_filename(filename, split_char=FILENAME_SPLIT,
                              mass_regex=MASS_REGEX):
    filename_elts = reversed(_filename_elts_list(filename, split_char))
    mass = _elt_from_felts(filename_elts, mass_regex)
    if mass is not None:
        return int(mass[1:])
    else:
        return None


def _filename_elts_list(filename, split_char):
    ext_index = filename.rfind('.')
    dir_index = filename.rfind('/')
    if dir_index != -1:
        filename_woext = filename[dir_index:ext_index]
    else:
        filename_woext = filename[:ext_index]
    return filename_woext.split(split_char)


def _elt_from_felts(felts, elt_regex):
    for elt in felts:
        m = re.match(elt_regex, elt)
        if m is not None and m.group(0) == elt:
            return elt
    else:
        return None


def _fname_without_extension(fname):
    index = fname.rfind('.')
    return fname[:index]


class SourcesDirDoesNotExistException(Exception):
    pass


def make_ans_file(file_path,
                  sp_file,
                  num_nucleons,
                  interaction_name='usdb',
                  option='lpe', neig=0,
                  restriction='n',
                  num_protons=8,
                  j_min=0.0, j_max=4.0, j_del=1.0,
                  parity=0,
                  end_option='st',
                  lines=LINES,
                  sep=SEP,
                  nl='\n'):
    """Create a .ans file with the given specifications
    :param nl:
    :param sep:
    :param lines:
    :param end_option:
    :param parity:
    :param j_del:
    :param j_max:
    :param j_min:
    :param num_protons:
    :param restriction:
    :param neig:
    :param option:
    :param interaction_name:
    :param num_nucleons:
    :param sp_file:
    :param file_path:
    """
    ans_str = nl.join(lines) % (sep,
                                option, neig,
                                sp_file,
                                restriction,
                                interaction_name,
                                num_protons,
                                num_nucleons,
                                j_min, j_max, j_del,
                                parity,
                                sep,
                                end_option)
    f = open(file_path, 'w')
    f.writelines(ans_str)
    f.close()


def main():
    do_all_calculations(
        arange=range(4, 20),
        zrange=[4, 5, 6],
        force=False
    )


main()
