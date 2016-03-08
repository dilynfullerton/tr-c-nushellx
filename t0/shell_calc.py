#!/usr/bin/python
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
DPATH_MAIN = getcwd()
DPATH_SOURCES = path.join(DPATH_MAIN, 'sources')
DPATH_RESULTS = path.join(DPATH_MAIN, 'results')
DNAME_USDB = 'usdb'
FNAME_MODEL_SPACE = 'imsrg'
FNAME_MODEL_SPACE_USDB = 'sd'

# file parsing
RGX_MASS_NUM = 'A\d+'
CHR_FILENAME_SPLIT = '_'
RGX_INT = '.*int'
RGX_ANS = '.*ans'
RGX_BAT = RGX_MASS_NUM + '\.bat'


def do_all_calculations(arange, zrange, **kwargs):
    for z in zrange:
        arange0 = list(filter(lambda a: a >= z, arange))
        make_results_dir(a_range=arange0, z=z, **kwargs)
        make_usdb_dir(a_range=arange0, z=z, **kwargs)
    return do_calculations(**kwargs)


def do_calculations(d=DPATH_MAIN,
                    dirpath_results=DPATH_RESULTS,
                    regex_ans=RGX_ANS,
                    regex_bat=RGX_BAT,
                    force=False):
    for root, dirs, files in walk(dirpath_results):
        ans = _get_file(files, regex_ans)
        bat = _get_file(files, regex_bat)
        if ans is not None:  # There is a *.ans file
            if bat is None or force:
                chdir(root)
                call(['shell', '%s' % ans])
    chdir(d)
    for root, dirs, files in walk(dirpath_results):
        bat = _get_file(files, regex_bat)
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


def files_with_ext_in_directory(dirpath, extension):
    """Returns a list of the filenames of all the files_INT in the given
    directory with the given extension
    :param extension:
    :param dirpath: """
    return list(glob.glob(path.join(dirpath, '*' + extension)))


def _get_file(list_of_fname, regex=RGX_ANS):
    for f in list_of_fname:
        if re.match(regex, f) is not None:
            return f
    else:
        return None


def make_usdb_dir(a_range, z,
                  d=DPATH_MAIN,
                  dirpath_results=DPATH_RESULTS,
                  dirname_usdb=DNAME_USDB,
                  fname_model_space=FNAME_MODEL_SPACE_USDB,
                  force=False):
    a_range = list(filter(lambda a: a >= z, a_range))
    results_subdir = path.join(dirpath_results, 'Z%d' % z)
    usdb_dirpath = path.join(results_subdir, dirname_usdb)
    if not path.exists(usdb_dirpath):
        mkdir(usdb_dirpath)
    for mass_num in a_range:
        dirname = path.join(usdb_dirpath, 'A%d' % mass_num)
        if not path.exists(dirname):
            mkdir(dirname)
        # link .sp file
        sp_filename = '%s.sp' % fname_model_space
        sp_file_path = path.join(dirname, sp_filename)
        if not path.exists(sp_file_path):
            link(path.join(d, sp_filename), sp_file_path)
        # create .ans file
        ans_filename = 'A%d.ans' % mass_num
        ans_file_path = path.join(dirname, ans_filename)
        if not path.exists(ans_file_path) or force is True:
            if mass_num % 2 == 0:  # even
                make_ans_file(file_path=ans_file_path,
                              sp_file=fname_model_space,
                              num_nucleons=mass_num,
                              interaction_name='usdb',
                              num_protons=z)
            else:
                make_ans_file(file_path=ans_file_path,
                              sp_file=fname_model_space,
                              num_nucleons=mass_num,
                              interaction_name='usdb',
                              num_protons=z,
                              j_min=0.5, j_max=3.5, j_del=1.0)


def make_results_dir(a_range, z,
                     d=DPATH_MAIN,
                     dirpath_sources=DPATH_SOURCES,
                     dirpath_results=DPATH_RESULTS,
                     fname_model_space=FNAME_MODEL_SPACE,
                     regex_int=RGX_INT,
                     force=False):
    """Copy all of the directories from the sources_dir into the results_dir
    recursively, but when encountering a *.int file, make a directory for the
    file (according to its name) and copy the file into the directory with a
    short name. Also, for each directory to which a *.int file is copied the
    given model space is linked and a *.ans file is generated.
    :param d: Main working directory
    :param dirpath_sources: Directory in which source files are stored
    :param dirpath_results: Directory in which results are to be generated
    :param fname_model_space: Name of the model space .sp file
    :param z: Z number for which to make results
    :param regex_int: Regular expression that matches interaction files
    :param force: If true, force making of .ans file, even if one
    already exists
    :param a_range: Mass range for which to create directories. If None,
    will do for all files.
    """
    a_range = list(filter(lambda a: a >= z, a_range))

    results_subdir = path.join(dirpath_results, 'Z%d' % z)
    if not path.exists(dirpath_sources):
        raise SourcesDirDoesNotExistException()
    if not path.exists(results_subdir):
        mkdir(results_subdir)

    todo_sources = deque([dirpath_sources])
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
        for ff in filter(lambda f: re.match(regex_int, f) is not None,
                         files):
            mass_num = mass_number_from_filename(ff)
            if a_range is not None and mass_num not in a_range:
                continue
            new_dir = path.join(cwd_results, _fname_without_extension(ff))
            if not path.exists(new_dir):
                mkdir(new_dir)
            # link .int file
            interaction_name = 'A%d' % mass_num
            interaction_file_path = path.join(new_dir,
                                              interaction_name + '.int')
            if not path.exists(interaction_file_path):
                link(path.join(cwd_sources, ff), interaction_file_path)
            # link .sd file
            sp_filename = '%s.sp' % fname_model_space
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
                                  sp_file=fname_model_space,
                                  restriction='n',
                                  interaction_name=interaction_name,
                                  num_protons=z,
                                  num_nucleons=mass_num,
                                  j_min=0.0, j_max=4.0, j_del=1.0,
                                  parity=0)
                else:
                    make_ans_file(file_path=ans_file_path,
                                  option='lpe', neig=0,
                                  sp_file=fname_model_space,
                                  restriction='n',
                                  interaction_name=interaction_name,
                                  num_protons=z,
                                  num_nucleons=mass_num,
                                  j_min=0.5, j_max=3.5, j_del=1.0,
                                  parity=0)


def mass_number_from_filename(filename,
                              split_char=CHR_FILENAME_SPLIT,
                              regex_mass_num=RGX_MASS_NUM):
    filename_elts = reversed(_filename_elts_list(filename, split_char))
    mass = _elt_from_felts(filename_elts, regex_mass_num)
    if mass is not None:
        return int(mass[1:])
    else:
        return None


def _filename_elts_list(fname, split_char):
    ext_index = fname.rfind('.')
    dir_index = fname.rfind('/')
    if dir_index != -1:
        filename_woext = fname[dir_index:ext_index]
    else:
        filename_woext = fname[:ext_index]
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
        arange=range(4, 7),
        zrange=[4, 5, 6],
        force=False
    )


if __name__ == "__main__":
    main()
