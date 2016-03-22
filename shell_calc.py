#!/usr/bin/python
"""shell_calc.py

To run as a script:

    $ shell_calc.py [-[Ff]] ncomponent formalism Amin [Amax] Zmin [Zmax]

For each valid (A, Z) pair in the range defined by [Amin, Amax] x [Zmin, Zmax],
where [.,.] signifies an inclusive interval over the positive integers,
runs the nushellx shell calculation for each matching interaction file in
SOURCE directory as well as usdb.

If -F or -f precedes the arguments, forces recalculation even if outfiles
    already exist in the RESULTS directory.
If 4 arguments are given,
    assumes these are ncomponent formalism Amin Zmin.
    Calculation range = {(Amin, Zmin)}, if Amin >= 2*Zmin and Amin, Zmin are
    positive integers.
If 5 arguments are given,
    assumes these are ncomponent formalism Amin Amax Zmin.
    Calculation range = the set of (A, Z) in [Amin, Amax] x [Zmin],
    where A >= 2*Z and A, Z are positive integers.
If 6 arguments are given,
    assumes these are ncomponent formalism Amin Amax Zmin Zmax.
    Calculation range = set of (A, Z) in [Amin, Amax] x [Zmin, Zmax],
    where A >= 2*Z and A, Z are positive integers.
"""
from __future__ import division

import glob
import re
from collections import deque
from math import floor
from os import getcwd, path, walk, mkdir, link, chdir, rmdir, listdir, remove
from subprocess import Popen, PIPE
from sys import argv, stdout
from threading import Thread, activeCount

# CONSTANTS
# .ans file
FLINES_FMT_ANS = [
    '--------------------------------------------------',
    '%-3s,   %-2d            ! option (lpe or lan), neig (zero=10)',
    '%-10s           ! model space (*.sp) name (a8)',
    '%s                    ! any restrictions (y/n)',
    '%-10s           ! interaction (*.int) name (a8)',
    '%3d                  ! number of protons',
    '%3d                  ! number of nucleons',
    '%4.1f,%4.1f,%4.1f,      ! min J, max J, del J',
    '%3d                  ! parity (0 for +) (1 for -) (2 for both)',
    '--------------------------------------------------',
    '%-3s                  ! option'
]

# directories
DPATH_MAIN = getcwd()
DPATH_SOURCES = path.join(DPATH_MAIN, 'sources')
DPATH_RESULTS = path.join(DPATH_MAIN, 'results')
DPATH_TEMPLATES = path.join(DPATH_MAIN, 'templates')
DNAME_FMT_Z = 'Z%d'
DNAME_USDB = 'usdb'

# file names
FNAME_MODEL_SPACE_P_N = 'tmp-p-n'
FNAME_MODEL_SPACE_P_PN = 'tmp-p-pn'
FNAME_MODEL_SPACE_SD_N = 'tmp-sd-n'
FNAME_MODEL_SPACE_SD_PN = 'tmp-sd-pn'
FNAME_MODEL_SPACE_OUT = 'aaa'
FNAME_MODEL_SPACE_USDB = 'sd'
_FNAME_SHELL_STDOUT = '__stdout_shell__.txt'
_FNAME_SHELL_STDERR = '__stderr_shell__.txt'
_FNAME_BAT_STDOUT = '__stdout_bat__.txt'
_FNAME_BAT_STDERR = '__stderr_bat__.txt'

# file parsing
RGX_MASS_NUM = 'A\d+'
_RGX_FNAME_INT = '.*\.int'
_RGX_FNAME_ANS = '.*\.ans'
_RGX_FNAME_BAT = RGX_MASS_NUM + '\.bat'
CHR_FILENAME_SPLIT = '_'

# printing
WIDTH_TERM = 79
WIDTH_PROGRESS_BAR = 48
STR_PROGRESS_BAR = '  Progress: %3d/%-3d '
STR_FMT_PROGRESS_HEAD = 'Doing shell calculation for Z = %d'

# shells
S_SHELL = frozenset(range(1, 4))
P_SHELL = frozenset(range(4, 16))
SD_SHELL = frozenset(range(16, 40))
MAP_SHELL_TO_MODEL_SPACES = {
    S_SHELL: (FNAME_MODEL_SPACE_P_N, FNAME_MODEL_SPACE_P_PN),
    # todo fix above value
    P_SHELL: (FNAME_MODEL_SPACE_P_N, FNAME_MODEL_SPACE_P_PN),
    SD_SHELL: (FNAME_MODEL_SPACE_SD_N, FNAME_MODEL_SPACE_SD_PN)
}

# threading
MAX_OPEN_THREADS = 12


class SourcesDirDoesNotExistException(Exception):
    pass


class NoAvailableModelSpaceException(Exception):
    pass


def _mass_number_from_filename(
        filename, split_char=CHR_FILENAME_SPLIT,
        regex_mass_num=RGX_MASS_NUM
):
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


def _make_ans_file(
        file_path, sp_file, num_nucleons, interaction_name='usdb',
        option='lpe', neig=0, restriction='n', num_protons=8,
        j_min=0.0, j_max=4.0, j_del=1.0, parity=0, end_option='st',
        lines=FLINES_FMT_ANS, nl='\n'
):
    """Create a .ans file with the given specifications
    :param nl: line separator
    :param lines: list of unformatted file lines
    :param end_option: final option to execute
    :param parity: parity
    :param j_del: difference between angular momentum values
    :param j_max: maximum j
    :param j_min: minimum j
    :param num_protons: proton number (z)
    :param restriction: restrictions?
    :param neig: something
    :param option: regular execution option
    :param interaction_name: name of interaction file
    :param num_nucleons: number of nucleons (A)
    :param sp_file: model space file name
    :param file_path: path to ans file
    """
    ans_str = nl.join(lines) % (
        option, neig, sp_file, restriction, interaction_name,
        num_protons, num_nucleons, j_min, j_max, j_del, parity, end_option)
    f = open(file_path, 'w')
    f.writelines(ans_str)
    f.close()


def _get_model_space(a, n_component=True,
                     shell_sp_map=MAP_SHELL_TO_MODEL_SPACES):
    """Returns the filename of the model space associated with the given
    mass number

    :param a: Mass number
    :param n_component: If True, assume a 2-component space (protons and
    neutrons), else assume a 1-component space (neutrons only)
    :param shell_sp_map: Map from shell to model space filename
    :return: Filename associated with the given A value
    :raises: NoAvailableModelSpaceException if there is no model space
    file for the given A value
    """
    for shell, sp in shell_sp_map.iteritems():
        if a in shell:
            return sp[n_component-1]
    else:
        if n_component == 2:
            s = ' pn'
        elif n_component == 1:
            s = ' n'
        else:
            s = ''
        raise NoAvailableModelSpaceException(
            'No%s model space available for A = %d' % (s, a)
        )


def _make_sp_file(src, dst, formalism):
    # read src file into list
    f = open(src, 'r')
    inlines = f.readlines()
    f.close()
    # apply replacements
    outlines = list()
    replacement_map = {'<<FORMALISM>>': formalism}
    for k, v in replacement_map.iteritems():
        for line in inlines:
            if k in line:
                outlines.append(line.replace(k, str(v)))
            else:
                outlines.append(line)
    # write new lines into dst
    f = open(dst, 'w')
    f.writelines(outlines)
    f.close()


def make_results_dir(
        a_range, z, ncomponent, formalism,
        dirpath_sources=DPATH_SOURCES,
        dirpath_results=DPATH_RESULTS,
        dirpath_templates=DPATH_TEMPLATES,
        _dname_fmt_z=DNAME_FMT_Z,
        _fname_model_space_out=FNAME_MODEL_SPACE_OUT,
        _regex_int=_RGX_FNAME_INT,
        force=False
):
    """Copy all of the directories from the sources_dir into the results_dir
    recursively, but when encountering a *.int file, make a directory for the
    file (according to its name) and copy the file into the directory with a
    short name. Also, for each directory to which a *.int file is copied the
    given model space is linked and a *.ans file is generated.
    :param a_range: Mass range for which to create directories. If None,
    will do for all files.
    :param z: Z number for which to make results
    :param ncomponent: Number of particle components in system.
        neutrons       -> 1
        proton-neutron -> 2
    :param formalism: formalism under which to run NuShellX.
        't'  -> isospin formalism
        'pn' -> proton-neutron formalism
    :param dirpath_sources: Directory in which source files are stored
    :param dirpath_results: Directory in which results are to be generated
    :param dirpath_templates: Directory in which *.sp templates files
    are housed
    :param _dname_fmt_z: Results subdirectory name template that takes the
    proton number (and integer) as its sole argument
    :param _fname_model_space_out: name of the model space file to write
    :param _regex_int: Regular expression that matches interaction files
    :param force: If true, force making of .ans file, even if one
    already exists
    """
    # only use A values greater than or equal to Z
    a_range = list(filter(lambda a: a >= z, a_range))

    results_subdir = path.join(dirpath_results, _dname_fmt_z % z)
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
        for ff in filter(lambda f: re.match(_regex_int, f) is not None, files):
            mass_num = _mass_number_from_filename(ff)
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
            elif force:
                remove(interaction_file_path)
                link(path.join(cwd_sources, ff), interaction_file_path)
            # write .sp file
            fname_model_space = _get_model_space(a=mass_num,
                                                 n_component=ncomponent)
            sp_filename = '%s.sp' % fname_model_space
            sp_filename_fin = '%s.sp' % _fname_model_space_out
            sp_path_src = path.join(dirpath_templates, sp_filename)
            sp_path_dst = path.join(new_dir, sp_filename_fin)
            if force or not path.exists(sp_path_dst):
                _make_sp_file(src=sp_path_src, dst=sp_path_dst,
                              formalism=formalism)
            # create .ans file
            ans_filename = 'A%d.ans' % mass_num
            ans_file_path = path.join(new_dir, ans_filename)
            if force or not path.exists(ans_file_path):
                if mass_num % 2 == 0:  # even
                    j_min, j_max, parity = 0.0, 4.0, 0
                else:  # odd
                    j_min, j_max, parity = 0.5, 3.5, 1
                _make_ans_file(
                    file_path=ans_file_path,
                    sp_file=_fname_model_space_out,
                    interaction_name=interaction_name,
                    option='lpe', neig=0, restriction='n',
                    num_protons=z, num_nucleons=mass_num,
                    j_min=j_min, j_max=j_max, j_del=1.0, parity=parity)


def make_usdb_dir(
        a_range, z,
        dirpath_results=DPATH_RESULTS,
        dirpath_templates=DPATH_TEMPLATES,
        _dirname_usdb=DNAME_USDB,
        _fname_model_space=FNAME_MODEL_SPACE_USDB,
        _usdb_arange=SD_SHELL,
        force=False
):
    a_range = list(filter(lambda a: a >= z, a_range))
    results_subdir = path.join(dirpath_results, 'Z%d' % z)
    usdb_dirpath = path.join(results_subdir, _dirname_usdb)
    if not path.exists(usdb_dirpath):
        mkdir(usdb_dirpath)
    for mass_num in a_range:
        if mass_num not in _usdb_arange:
            continue
        dirname = path.join(usdb_dirpath, 'A%d' % mass_num)
        if not path.exists(dirname):
            mkdir(dirname)
        # link .sp file
        sp_filename = '%s.sp' % _fname_model_space
        sp_file_path = path.join(dirname, sp_filename)
        if not path.exists(sp_file_path):
            link(path.join(dirpath_templates, sp_filename), sp_file_path)
        elif force:
            remove(sp_file_path)
            link(path.join(dirpath_templates, sp_filename), sp_file_path)
        # create .ans file
        ans_filename = 'A%d.ans' % mass_num
        ans_file_path = path.join(dirname, ans_filename)
        if force or not path.exists(ans_file_path):
            if mass_num % 2 == 0:  # even
                j_min, j_max, parity = 0.0, 4.0, 0
            else:
                j_min, j_max, parity = 0.5, 3.5, 0
            _make_ans_file(
                file_path=ans_file_path, sp_file=_fname_model_space,
                num_nucleons=mass_num, num_protons=z,
                interaction_name='usdb',
                j_min=j_min, j_max=j_max, j_del=1.0, parity=parity, )


def remove_empty_directories(root, remove_root=False):
    item_names = listdir(root)
    item_paths = [path.join(root, item) for item in item_names]
    subdir_paths = list(filter(lambda p: path.isdir(p), item_paths))
    for sd in subdir_paths:
        remove_empty_directories(root=sd, remove_root=True)
    item_names = listdir(root)
    if len(item_names) == 0 and remove_root:
        rmdir(root)


def _get_file(list_of_fname, regex=_RGX_FNAME_ANS):
    for f in list_of_fname:
        if re.match(regex, f) is not None:
            return f
    else:
        return None


def _files_with_ext_in_directory(dirpath, extension):
    """Returns a list of the filenames of all the files_INT in the given
    directory with the given extension
    :param extension: file extension
    :param dirpath: path to directory
    """
    return list(glob.glob(path.join(dirpath, '*' + extension)))


def _calc_has_been_done(dirpath):
    return len(_files_with_ext_in_directory(dirpath, '.lpt')) > 1


def _shell_calculation(
        root, fname_ans, verbose,
        _fname_stdout=_FNAME_SHELL_STDOUT, _fname_stderr=_FNAME_SHELL_STDERR
):
    main_dir = getcwd()
    chdir(root)
    args = ['shell', '%s' % fname_ans]
    if not verbose:
        p = Popen(args=args, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        fout = open(path.join(root, _fname_stdout), 'w')
        fout.write(out)
        fout.close()
        ferr = open(path.join(root, _fname_stderr), 'w')
        ferr.write(err)
        ferr.close()
    else:
        p = Popen(args=args)
        p.wait()
    chdir(main_dir)


def _do_shell_calculation(
        files, root, force, verbose,
        _rgx_ans=_RGX_FNAME_ANS, _rgx_bat=_RGX_FNAME_BAT
):
    fname_ans = _get_file(files, _rgx_ans)
    fname_bat = _get_file(files, _rgx_bat)
    if fname_ans is not None:  # There is a *.ans file
        if fname_bat is None or force:
            _shell_calculation(root=root, fname_ans=fname_ans, verbose=verbose)


def _bat_calculation(
        root, fname_bat, verbose,
        _fname_stdout=_FNAME_BAT_STDOUT, _fname_stderr=_FNAME_BAT_STDERR
):
    main_dir = getcwd()
    chdir(root)
    args = ['source', '%s' % fname_bat]
    if not verbose:
        try:
            p = Popen(args=args, stdout=PIPE, stderr=PIPE)
        except OSError:
            p = Popen(args=' '.join(args), shell=True,
                      stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        fout = open(path.join(root, _fname_stdout), 'w')
        fout.write(out)
        fout.close()
        ferr = open(path.join(root, _fname_stderr), 'w')
        ferr.write(err)
        ferr.close()
    else:
        try:
            p = Popen(args=args)
        except OSError:
            p = Popen(args=' '.join(args), shell=True)
        p.wait()
    chdir(main_dir)


def _do_bat_calculation(
        root, files, force, verbose, _rgx_bat=_RGX_FNAME_BAT
):
    fname_bat = _get_file(files, _rgx_bat)
    if fname_bat is not None:
        if not _calc_has_been_done(root) or force is True:
            _bat_calculation(root=root, fname_bat=fname_bat, verbose=verbose)


def _print_progress(
        completed, total, end=False,
        bar_len=WIDTH_PROGRESS_BAR,
        total_width=WIDTH_TERM,
        text_fmt=STR_PROGRESS_BAR
):
    if total > 0:
        text = text_fmt % (floor(completed), total)
        p = completed / total
        bar_fill = int(floor(p * bar_len))
        progress_bar = '[' + '#'*bar_fill + ' '*(bar_len - bar_fill) + ']'
        sp_fill_len = total_width - len(text) - len(progress_bar)
        if sp_fill_len < 0:
            sp_fill_len = 0
        line = '\r' + text + ' '*sp_fill_len + progress_bar
        if end:
            line += '\n'
        stdout.write(line)
        stdout.flush()


def _shell_and_bat(root, files, force, verbose):
    _do_shell_calculation(root=root, files=files,
                          force=force, verbose=verbose)
    new_files = listdir(root)
    _do_bat_calculation(root=root, files=new_files,
                        force=force, verbose=verbose)


def _do_calculation_t(
        todo_walk, z, force, progress,
        _max_open_threads=MAX_OPEN_THREADS,
        _str_fmt_prog=STR_FMT_PROGRESS_HEAD
):
    def _r(root_, files_):
        _shell_and_bat(root=root_, files=files_, force=force, verbose=False)
    threads_opened = deque()
    todo_list = list(todo_walk)
    jobs_completed = 0
    jobs_total = len(todo_list)
    if progress:
        print _str_fmt_prog % z
    while len(todo_list) > 0 or len(threads_opened) > 0:
        if progress:
            _print_progress(jobs_completed, jobs_total)
        # if room, start new threads
        while activeCount() < _max_open_threads and len(todo_list) > 0:
            root, files = todo_list.pop()
            t = Thread(target=_r, args=(root, files))
            t.start()
            threads_opened.append(t)
        # wait for completion of first thread
        t = threads_opened.popleft()
        t.join()
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def _do_calculation(todo_walk, z, force, verbose, progress,
                    _str_fmt_prog=STR_FMT_PROGRESS_HEAD):
    jobs_completed = 0
    jobs_total = len(todo_walk)
    if progress:
        print _str_fmt_prog % z
    for root, files in todo_walk:
        # do shell calculation
        if progress:
            _print_progress(jobs_completed, jobs_total)
        _shell_and_bat(root=root, files=files, force=force, verbose=verbose)
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def do_calculations(
        a_range, z,
        dirpath_main=DPATH_MAIN,
        dirpath_results=DPATH_RESULTS,
        dirname_fmt_z=DNAME_FMT_Z,
        _rgx_bat=_RGX_FNAME_BAT, _rgx_int=_RGX_FNAME_INT,
        force=False, verbose=False, progress=True, threading=True
):
    dirpath_z = path.join(dirpath_results, dirname_fmt_z % z)
    todo = list()
    for root, dirs, files in walk(dirpath_z):
        # if the mass number is not in the range, do not do calculation
        fname_int = _get_file(files, _rgx_int)
        fname_bat = _get_file(files, _rgx_bat)
        if fname_int is not None and (force or fname_bat is None):
            a = _mass_number_from_filename(filename=fname_int)
            if a in a_range:
                todo.append((root, files))
    if threading and len(todo) > 1:
        _do_calculation_t(todo_walk=todo, z=z, force=force, progress=progress)
    else:
        _do_calculation(
            todo_walk=todo, z=z,
            force=force, progress=progress, verbose=verbose,
        )
    chdir(dirpath_main)
    return 1


def do_all_calculations(
        arange, zrange, n_component=2, formalism='pn',
        dirpath_results=DPATH_RESULTS, **kwargs
):
    zrange = list(filter(lambda z0: z0 >= 1, zrange))
    for z in zrange:
        arange0 = list(filter(lambda a: a >= z, arange))
        make_results_dir(a_range=arange0, z=z,
                         ncomponent=n_component,
                         formalism=formalism, **kwargs)
        make_usdb_dir(a_range=arange0, z=z, **kwargs)
        remove_empty_directories(dirpath_results, remove_root=False)
        do_calculations(a_range=arange0, z=z, **kwargs)


if __name__ == "__main__":
    user_args = argv[1:]
    if re.match('^-f$', user_args[0].lower()):
        force0 = True
        user_args = user_args[1:]
    else:
        force0 = False
    if len(user_args) == 4:
        ncomponent0 = int(user_args[0])
        formalism0 = user_args[1]
        amin, zmin = [int(x) for x in user_args[2:]]
        do_all_calculations(arange=[amin], zrange=[zmin],
                            n_component=ncomponent0, formalism=formalism0,
                            force=force0)
    elif len(user_args) == 5:
        ncomponent0 = int(user_args[0])
        formalism0 = user_args[1]
        amin, amax, zmin = [int(x) for x in user_args[2:]]
        do_all_calculations(arange=range(amin, amax+1), zrange=[zmin],
                            n_component=ncomponent0, formalism=formalism0,
                            force=force0)
    elif len(user_args) == 6:
        ncomponent0 = int(user_args[0])
        formalism0 = user_args[1]
        amin, amax, zmin, zmax = [int(x) for x in user_args[2:]]
        do_all_calculations(arange=range(amin, amax+1),
                            zrange=range(zmin, zmax+1),
                            n_component=ncomponent0, formalism=formalism0,
                            force=force0)
    else:
        print ('User entered %d arguments. ' % (len(user_args),) +
               'shell_calc.py requires 4-6 arguments.')
