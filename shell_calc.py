#!/usr/bin/python
"""shell_calc.py

To run as a script:

    $ shell_calc.py [-[Ff]] nshell ncomponent formalism Amin [Amax] Zmin [Zmax]

For each valid (A, Z) pair in the range defined by [Amin, Amax] x [Zmin, Zmax],
where [.,.] signifies an inclusive interval over the positive integers,
runs the nushellx shell calculation for each matching interaction file in
SOURCE directory as well as usdb.

If -F or -f precedes the arguments, forces recalculation even if outfiles
    already exist in the RESULTS directory.
If 5 arguments are given,
    assumes these are nshell ncomponent formalism Amin Zmin.
    Calculation range = {(Amin, Zmin)}, if Amin >= 2*Zmin and Amin, Zmin are
    positive integers.
If 6 arguments are given,
    assumes these are nshell ncomponent formalism Amin Amax Zmin.
    Calculation range = the set of (A, Z) in [Amin, Amax] x [Zmin],
    where A >= 2*Z and A, Z are positive integers.
If 7 arguments are given,
    assumes these are nshell ncomponent formalism Amin Amax Zmin Zmax.
    Calculation range = set of (A, Z) in [Amin, Amax] x [Zmin, Zmax],
    where A >= 2*Z and A, Z are positive integers.
"""
from __future__ import division

import glob
import re
from Queue import Queue
from collections import deque
from math import floor
from os import getcwd, path, walk, mkdir, link, rmdir, listdir, remove, makedirs
from subprocess import Popen, PIPE
from sys import argv, stdout
from threading import Thread, currentThread

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
JMIN_EVEN = 0.0
JMAX_EVEN = 4.0
PARITY_EVEN = 0
JMIN_ODD = JMIN_EVEN + 0.5
JMAX_ODD = JMAX_EVEN - 0.5
PARITY_ODD = 1

# directories
DPATH_MAIN = getcwd()
DPATH_SOURCES = path.join(DPATH_MAIN, 'sources')
DPATH_RESULTS = path.join(DPATH_MAIN, 'results')
DPATH_TEMPLATES = path.join(DPATH_MAIN, 'templates')
DNAME_FMT_Z = 'Z%d'
DNAME_USDB = 'usdb'

# file names
FNAME_MODEL_SPACE_P_N = 'pn'
FNAME_MODEL_SPACE_P_PN = 'ppn'
FNAME_MODEL_SPACE_SD_N = 'sdn'
FNAME_MODEL_SPACE_SD_PN = 'sdpn'
FNAME_MODEL_SPACE_USDB = 'sd'
_FNAME_STDOUT_SHELL = '__stdout_shell__.txt'
_FNAME_STDERR_SHELL = '__stderr_shell__.txt'
_FNAME_STDOUT_BAT = '__stdout_bat__.txt'
_FNAME_STDERR_BAT = '__stderr_bat__.txt'

# file parsing
RGX_MASS_NUM = 'A\d+'
_RGX_FNAME_INT = '(A\d+|usdb)\.int'
_RGX_FNAME_ANS = RGX_MASS_NUM + '\.ans'
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
MAX_OPEN_THREADS = 10


class SourcesDirDoesNotExistException(Exception):
    pass


class NoAvailableModelSpaceException(Exception):
    pass


def _mass_number_from_filename(
        filename, split_char=CHR_FILENAME_SPLIT, regex_mass_num=RGX_MASS_NUM):
    """Given a file name, which is formatted with "elements" or "blocks"
    separated by a split_char, retrieves the one that specifies the mass number
    and returns the mass number.
    :param filename: file name of a NuShellX *.int file
    :param split_char: character that separates filename "elements". Default is
    '_' (underscore)
    :param regex_mass_num: regular expression that matches the mass number
    in the filename. Example: "A\d+"
    :return: mass number (an integer)
    """
    filename_elts = reversed(_filename_elts_list(filename, split_char))
    mass = _elt_from_felts(filename_elts, regex_mass_num)
    if mass:
        return int(mass[1:])
    else:
        return None


def _filename_elts_list(fname, split_char):
    """Returns a list of filename "elements," which are the pieces of the
    filename separated by spit_char.
    Extension of the filename is removed.
    If the filename is actually a file path, directories are removed leaving
    only the filename
    :param fname: name of the file
    :param split_char: character that splits the filename "elements" or
    "blocks". For example, '_' (underscore)
    :return: list of filename elements
    """
    ext_index = fname.rfind('.')
    dir_index = fname.rfind('/')
    if dir_index != -1:
        filename_without_ext = fname[dir_index:ext_index]
    else:
        filename_without_ext = fname[:ext_index]
    return filename_without_ext.split(split_char)


def _elt_from_felts(felts, elt_regex):
    """Returns the first of the elements in felts that matches elt_regex
    completely
    :param felts: list of filename elements (strings)
    :param elt_regex: regular expression that matches the desired element
    :return: first element in felts that matches elt_regex or None if a
    matching element is not found
    """
    for elt in felts:
        m = re.match(elt_regex, elt)
        if m and m.group(0) == elt:
            return elt
    else:
        return None


def _fname_without_extension(fname):
    """Given a filename with an extension, return the filename without the
    extension
    :param fname: a filename
    :return: given filename without everything following (and including) the
    rightmost period, if it exists
    """
    index = fname.rfind('.')
    if index != -1:
        return fname[:index]
    else:
        return fname


def _make_ans_file(
        file_path, sp_file, num_nucleons, interaction_name='usdb',
        option='lpe', neig=0, restriction='n', num_protons=8,
        j_min=0.0, j_max=4.0, j_del=1.0, parity=0, end_option='st',
        lines=FLINES_FMT_ANS,
):
    """Create a .ans file with the given specifications
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
    ans_str = '\n'.join(lines) % (
        option, neig, sp_file, restriction, interaction_name,
        num_protons, num_nucleons, j_min, j_max, j_del, parity, end_option
    )
    f = open(file_path, 'w')
    f.writelines(ans_str)
    f.close()


def _get_model_space(
        a, n_component=True, shell_sp_map=MAP_SHELL_TO_MODEL_SPACES):
    """Returns the filename of the model space associated with the given
    mass number
    :param a: mass number (A)
    :param n_component: (1 -> neutrons, 2 -> protons and neutrons)
    :param shell_sp_map: Map from shell to model space filename
    :return: filename associated with the given A value
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
            'No%s model space available for A = %d' % (s, a))


def _make_sp_file(src, dst, formalism):
    """Using the template file given by src, write a model space *.sp file at
    dst with the given formalism
    :param src: path to the appropriate template *.sp file
    :param dst: path to the destination *.sp file
    :param formalism: t -> isopsin, pn -> proton/neutron
    """
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


class UnknownParityException(Exception):
    pass


# todo make me general
def _get_parity(a, nshell):
    """Given a mass number and shell, return the parity
    :param a: mass number (A)
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :raises UnknownParityException: raised if the parity is unknown for the
    given shell
    """
    if nshell == 1:
        return a % 2
    elif nshell == 2:
        return 0
    else:
        raise UnknownParityException(
            '\nMethod of determining parity is not known for '
            'Nshell = %d' % nshell
        )


def make_results_dir(
        a_range, z, nshell, ncomponent, formalism,
        jmin_even=JMIN_EVEN, jmax_even=JMAX_EVEN,
        jmin_odd=JMIN_ODD, jmax_odd=JMAX_ODD,
        dpath_sources=DPATH_SOURCES,
        dpath_results=DPATH_RESULTS,
        dpath_templates=DPATH_TEMPLATES,
        _dname_fmt_z=DNAME_FMT_Z,
        _regex_int=_RGX_FNAME_INT,
        force=False
):
    """Copy all of the directories from the sources_dir into the results_dir
    recursively, but when encountering a *.int file, make a directory for the
    file (according to its name) and copy the file into the directory with a
    short name. Also, for each directory to which a *.int file is copied the
    given model space is linked and a *.ans file is generated.
    :param a_range: mass range for which to create directories. If None,
    will do for all files.
    :param z: Z number for which to make results
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param ncomponent: number of particle components in system.
        neutrons       -> 1
        proton-neutron -> 2
    :param formalism: formalism under which to run NuShellX.
        't'  -> isospin formalism
        'pn' -> proton-neutron formalism
    :param jmin_even: *.ans parameter jmin for even A
    :param jmax_even: *.ans parameter jmax for even A
    :param jmin_odd: *.ans parameter jmin for odd A
    :param jmax_odd: *.ans parameter jmax for odd A
    :param dpath_sources: directory in which source files are stored
    :param dpath_results: directory in which results are to be generated
    :param dpath_templates: directory in which *.sp templates files
    are housed
    :param _dname_fmt_z: results subdirectory name template that takes the
    proton number (and integer) as its sole argument
    :param _regex_int: regular expression that matches interaction files
    :param force: If true, force making of .ans file, even if one
    already exists
    """
    # only use A values greater than or equal to Z
    a_range = list(filter(lambda a: a >= z, a_range))
    results_subdir = path.join(dpath_results, _dname_fmt_z % z)
    if not path.exists(dpath_sources):
        raise SourcesDirDoesNotExistException()
    if not path.exists(results_subdir):
        makedirs(results_subdir)
    todo_sources = deque([dpath_sources])
    todo_results = deque([results_subdir])
    while len(todo_sources) > 0:
        cwd_sources = todo_sources.popleft()
        cwd_results = todo_results.popleft()
        root, dirs, files = walk(cwd_sources).next()
        for dname in dirs:
            next_sources = path.join(cwd_sources, dname)
            next_results = path.join(cwd_results, dname)
            if not path.exists(next_results):
                mkdir(next_results)
            todo_sources.append(next_sources)
            todo_results.append(next_results)
        for fname in filter(lambda f: re.match(_regex_int, f), files):
            mass_num = _mass_number_from_filename(fname)
            if a_range is not None and mass_num not in a_range:
                continue
            new_dir = path.join(cwd_results, _fname_without_extension(fname))
            if not path.exists(new_dir):
                mkdir(new_dir)
            # link .int file
            interaction_name = 'A%d' % mass_num
            fpath_int_dst = path.join(
                new_dir, interaction_name + '.int')
            if force or not path.exists(fpath_int_dst):
                if path.exists(fpath_int_dst):
                    remove(fpath_int_dst)
                link(path.join(cwd_sources, fname), fpath_int_dst)
            # write .sp file
            fname_model_space = _get_model_space(
                a=mass_num, n_component=ncomponent)
            fname_sp_src = '%s.sp' % fname_model_space
            fname_sp_dst = '%s.sp' % fname_model_space
            fpath_sp_src = path.join(dpath_templates, fname_sp_src)
            fpath_sp_dst = path.join(new_dir, fname_sp_dst)
            if force or not path.exists(fpath_sp_dst):
                _make_sp_file(
                    src=fpath_sp_src, dst=fpath_sp_dst, formalism=formalism)
            # create .ans file
            fpath_ans = path.join(new_dir, 'A%d.ans' % mass_num)
            if force or not path.exists(fpath_ans):
                if mass_num % 2 == 0:  # even
                    j_min, j_max = jmin_even, jmax_even
                    parity = _get_parity(a=mass_num, nshell=nshell)
                else:  # odd
                    j_min, j_max = jmin_odd, jmax_odd
                    parity = _get_parity(a=mass_num, nshell=nshell)
                _make_ans_file(
                    file_path=fpath_ans,
                    sp_file=fname_model_space,
                    interaction_name=interaction_name,
                    option='lpe', neig=0, restriction='n',
                    num_protons=z, num_nucleons=mass_num,
                    j_min=j_min, j_max=j_max, j_del=1.0, parity=parity
                )


def make_usdb_dir(
        a_range, z, nshell,
        jmin_even=JMIN_EVEN, jmax_even=JMAX_EVEN,
        jmin_odd=JMIN_ODD, jmax_odd=JMAX_ODD,
        dpath_results=DPATH_RESULTS,
        dpath_templates=DPATH_TEMPLATES,
        _dname_usdb=DNAME_USDB,
        _fname_model_space=FNAME_MODEL_SPACE_USDB,
        _usdb_arange=SD_SHELL,
        force=False
):
    """Make USDB directories for the given parameters. This is largely similar
    to make_results_dir, but in this case the interaction is always usdb.int
    :param a_range: mass range for which to create directories. If None,
    will do for all files.
    :param z: Z number for which to make results
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param jmin_even: *.ans parameter jmin for even A
    :param jmax_even: *.ans parameter jmax for even A
    :param jmin_odd: *.ans parameter jmin for odd A
    :param jmax_odd: *.ans parameter jmax for odd A
    :param dpath_results: directory in which results are to be generated
    :param dpath_templates: directory in which *.sp templates files
    are housed
    :param _dname_usdb: name of the usdb subdirectory
    :param _fname_model_space: name of the usdb model space file
    :param _usdb_arange: range of A values for which USDB can be done
    :param force: If true, force making of .ans file, even if one
    already exists
    :return:
    """
    a_range = list(filter(lambda a: a >= z, a_range))
    results_subdir = path.join(dpath_results, 'Z%d' % z)
    usdb_dirpath = path.join(results_subdir, _dname_usdb)
    if not path.exists(usdb_dirpath):
        makedirs(usdb_dirpath)
    for mass_num in a_range:
        if mass_num not in _usdb_arange:
            continue
        dname = path.join(usdb_dirpath, 'A%d' % mass_num)
        if not path.exists(dname):
            mkdir(dname)
        # link .sp file
        fname_sp = '%s.sp' % _fname_model_space
        fpath_sp = path.join(dname, fname_sp)
        if force or not path.exists(fpath_sp):
            if path.exists(fpath_sp):
                remove(fpath_sp)
            link(path.join(dpath_templates, fname_sp), fpath_sp)
        # create .ans file
        fpath_ans = path.join(dname, 'A%d.ans' % mass_num)
        if force or not path.exists(fpath_ans):
            if mass_num % 2 == 0:  # even
                j_min, j_max = jmin_even, jmax_even
                parity = _get_parity(a=mass_num, nshell=nshell)
            else:  # odd
                j_min, j_max = jmin_odd, jmax_odd
                parity = _get_parity(a=mass_num, nshell=nshell)
            _make_ans_file(
                file_path=fpath_ans, sp_file=_fname_model_space,
                num_nucleons=mass_num, num_protons=z, interaction_name='usdb',
                j_min=j_min, j_max=j_max, j_del=1.0, parity=parity,
            )


def remove_empty_directories(root, remove_root=False):
    """Starting from a given root, recursively remove any subdirectory that
    contains no files or directories (besides other empty directories)
    :param root: starting directory
    :param remove_root: if true, will remove the root if it is empty after
    all empty subdirectories are removed
    """
    item_names = listdir(root)
    item_paths = [path.join(root, item) for item in item_names]
    subdir_paths = list(filter(lambda p: path.isdir(p), item_paths))
    for sd in subdir_paths:
        remove_empty_directories(root=sd, remove_root=True)
    item_names = listdir(root)
    if len(item_names) == 0 and remove_root:
        rmdir(root)


# todo does this really need its own function?
def _get_file(list_of_fname, regex=_RGX_FNAME_ANS):
    """In a list of filenames, returns the file that matches the given
    regular expression.
    :param list_of_fname: list of filenames to search
    :param regex: regular expression
    :return: first file that matches regex, or None if no matching file is
    found
    """
    for f in list_of_fname:
        if re.match(regex, f):
            return f
    else:
        return None


def _files_with_ext_in_directory(dirpath, extension):
    """Returns a list of the filenames of all the files in the given
    directory with the given extension
    :param extension: file extension
    :param dirpath: path to directory
    """
    return list(glob.glob(path.join(dirpath, '*' + extension)))


def _calc_has_been_done(dpath):
    """Determine whether the shell calculation has completed successfully.
    If it has, there should be multiple files with the .lpt extension.
    :param dpath: path to the directory to search
    :return: true if there is more that one file with the extension '.lpt',
    false otherwise
    """
    return len(_files_with_ext_in_directory(dpath, '.lpt')) > 1


def _shell_calculation(
        root, fname_ans, verbose,
        _fname_stdout=_FNAME_STDOUT_SHELL, _fname_stderr=_FNAME_STDERR_SHELL
):
    """Run shell in root
    :param root: path to the directory in which to run shell
    :param fname_ans: name of the *.ans file specifying how the calculation is
    to be done
    :param verbose: if true, regular output of shell is outputted to stdout
    :param _fname_stdout: filename in which to save standard output of shell
    if verbose is false
    :param _fname_stderr: filename in which to save error output of shell
    if verbose is false
    :param return: returncode
    """
    args = ['shell', '%s' % fname_ans]
    if not verbose:
        p = Popen(args=args, stdout=PIPE, stderr=PIPE, cwd=root)
        out, err = p.communicate()
        if len(out) > 0:
            fout = open(path.join(root, _fname_stdout), 'w')
            fout.write(out)
            fout.close()
        if len(err) > 0:
            ferr = open(path.join(root, _fname_stderr), 'w')
            ferr.write(err)
            ferr.close()
    else:
        p = Popen(args=args, cwd=root)
        p.wait()
    return p.poll()


def _do_shell_calculation(
        files, root, force, verbose,
        _rgx_ans=_RGX_FNAME_ANS, _rgx_bat=_RGX_FNAME_BAT
):
    """Do shell calculation in root IF it has not already been done or force
    is True
    :param files: files in root
    :param root: directory in which to run shell
    :param force: if true, runs shell even if output files exist
    :param verbose: if true, regular output of shell is printed to stdout
    :param _rgx_ans: regular expression matching *.ans file
    :param _rgx_bat: regular expression matching *.bat file
    :return: returncode of shell calculation
    """
    fname_ans = _get_file(files, _rgx_ans)
    fname_bat = _get_file(files, _rgx_bat)
    if fname_ans is not None:  # There is a *.ans file
        if fname_bat is None or force:
            return _shell_calculation(
                root=root, fname_ans=fname_ans, verbose=verbose)


def _bat_calculation(
        root, fname_bat, verbose,
        _fname_stdout=_FNAME_STDOUT_BAT, _fname_stderr=_FNAME_STDERR_BAT
):
    """Source the given *.bat file in root
    :param root: directory in which to run calculation
    :param fname_bat: name of the *.bat file
    :param verbose: if true, prints standard output to stdout
    :param _fname_stdout: filename to which to write standard output if
    verbose is false
    :param _fname_stderr: filenaem to which to write error output if
    verbose is false
    :return: returncode
    """
    args = ['source', '%s' % fname_bat]
    if not verbose:
        try:
            p = Popen(args=args, stdout=PIPE, stderr=PIPE, cwd=root)
        except OSError:
            p = Popen(args=' '.join(args), shell=True,
                      stdout=PIPE, stderr=PIPE, cwd=root)
        out, err = p.communicate()
        if len(out) > 0:
            fout = open(path.join(root, _fname_stdout), 'w')
            fout.write(out)
            fout.close()
        if len(err) > 0:
            ferr = open(path.join(root, _fname_stderr), 'w')
            ferr.write(err)
            ferr.close()
    else:
        try:
            p = Popen(args=args, cwd=root)
        except OSError:
            p = Popen(args=' '.join(args), shell=True, cwd=root)
        p.wait()
    return p.poll()


def _do_bat_calculation(
        root, files, force, verbose, _rgx_bat=_RGX_FNAME_BAT
):
    """Do bat calculation IF it has not already been done or if the user
     forces
    :param root: directory in which to do calculation
    :param files: files in the directory
    :param force: if true, runs *.bat even if output files already exist
    :param verbose: if true, writes regular output to standard out
    :param _rgx_bat: regular expression that matches the *.bat file to run
    :return: returncode
    """
    fname_bat = _get_file(files, _rgx_bat)
    if fname_bat is not None:
        if force or not _calc_has_been_done(root):
            return _bat_calculation(
                root=root, fname_bat=fname_bat, verbose=verbose)


def _print_progress(
        completed, total, end=False,
        bar_len=WIDTH_PROGRESS_BAR,
        total_width=WIDTH_TERM,
        text_fmt=STR_PROGRESS_BAR
):
    """Writes a progress bar to stdout based on completed and total
    :param completed: jobs completed
    :param total: total jobs todo
    :param end: if true, adds a newline '\n' character to the end of output to
    allow for regular printing
    :param bar_len: width of the error bar
    :param total_width: total allowed width of output
    :param text_fmt: text to be printed beside the bar, which takes completed
    and total as format arguments
    """
    if not total > 0:
        return None
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
    """Do the shell and bat calculations
    :param root: directory in which to do the calculations
    :param files: names of files in the directory
    :param force: if true, forces redoing of calculations even if output
    files already exist
    :param verbose: if true, prints regular output to stdout; else, this is
    suppressed and written to files instead
    :return: exit codes for shell and bat calculations
    """
    ret1 = _do_shell_calculation(
        root=root, files=files, force=force, verbose=verbose)
    new_files = listdir(root)
    ret2 = _do_bat_calculation(
        root=root, files=new_files, force=force, verbose=verbose)
    return ret1, ret2


def _do_calculation_t(
        todo_walk, z, force, progress,
        _max_open_threads=MAX_OPEN_THREADS,
        _str_fmt_prog=STR_FMT_PROGRESS_HEAD
):
    """Run the calculations in concurrency with a maximum number of open
    threads.
    :param todo_walk: list of all pairs (root, files) for which the
    shell and bat calculations are to be done
    :param z: proton number (Z)
    :param force: if true, forces redoing shell and bat even if output files
    already exist
    :param progress: if true, shows a progress bar, indicating the number
    of closed threads
    :param _max_open_threads: maximum number of external threads to open
    :param _str_fmt_prog: string to display above the progress bar
    :return: true, if all jobs were completed; false otherwise
    """
    def _r(root_, files_, q):
        _shell_and_bat(root=root_, files=files_, force=force, verbose=False)
        q.put(currentThread())
    todo_list = list(todo_walk)
    active_list = list()
    done_queue = Queue()
    jobs_completed = 0
    jobs_total = len(todo_list)
    if progress and jobs_total > 0:
        print _str_fmt_prog % z
        _print_progress(jobs_completed, jobs_total)
    while len(todo_list) > 0 or len(active_list) > 0:
        # if room, start new threads
        while len(todo_list) > 0 and len(active_list) < _max_open_threads:
            root, files = todo_list.pop()
            t = Thread(target=_r, args=(root, files, done_queue))
            active_list.append(t)
            t.start()
        # remove any threads that have finished
        if not done_queue.empty():
            while not done_queue.empty():
                t = done_queue.get()
                t.join()
                active_list.remove(t)
                jobs_completed += 1
            if progress:
                _print_progress(jobs_completed, jobs_total)
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    return jobs_completed == jobs_total


def _do_calculation(todo_walk, z, force, verbose, progress,
                    _str_fmt_prog=STR_FMT_PROGRESS_HEAD):
    """Run shell and bat calculations sequentially
    :param todo_walk: list of (root, files) for which the calculations are to
    be done
    :param z: proton number (Z)
    :param force: if true, forces redoing calculations even if output files
    exist
    :param verbose: if true, writes output of shell and bat to stdout; else
    output is suppressed and written to files instead
    :param progress: if true, displays a progress bar (not compatible with
    verbose option)
    :param _str_fmt_prog: string to display above the progress bar
    :return: true if all jobs were completed; false otherwise
    """
    jobs_completed = 0
    jobs_total = len(todo_walk)
    progress = progress and not verbose
    if progress and jobs_total > 0:
        print _str_fmt_prog % z
    for root, files in todo_walk:
        # do shell calculation
        if progress:
            _print_progress(jobs_completed, jobs_total)
        _shell_and_bat(root=root, files=files, force=force, verbose=verbose)
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    return jobs_completed == jobs_total


def do_calculations(
        a_range, z,
        dpath_results=DPATH_RESULTS, dname_fmt_z=DNAME_FMT_Z,
        _rgx_bat=_RGX_FNAME_BAT, _rgx_int=_RGX_FNAME_INT,
        force=False, verbose=False, progress=True, threading=True,
):
    """For a given proton number, do all shell and bat calculations that have
    not been done in the results directory
    :param a_range: list of mass numbers for which to do calculations
    :param z: proton number (Z) for which to do calculations
    :param dpath_results: results directory
    :param dname_fmt_z: results subdirectory for the the given proton number,
    (to be formatted with that number)
    :param _rgx_bat: regular expression that matches the *.bat files
    :param _rgx_int: regular expression that matches the *.int files
    :param force: if true, redoes calculations even if output files already
    exist
    :param verbose: if true, write standard output to stdout; else this is
    suppressed and written to files
    :param progress: if true, prints a progress bar to stdout (not compatible
    with verbose)
    :param threading: if true, uses multithreading to speed up calculation
    :return: true if all calculations were completed successfully; false
    otherwise
    """
    progress = progress and not verbose
    dirpath_z = path.join(dpath_results, dname_fmt_z % z)
    todo = list()
    for root, dirs, files in walk(dirpath_z):
        # if the mass number is not in the range, do not do calculation
        fname_int = _get_file(files, _rgx_int)
        fname_bat = _get_file(files, _rgx_bat)
        if fname_int and (force or not fname_bat):
            a = _mass_number_from_filename(filename=fname_int)
            if a in a_range:
                todo.append((root, files))
    if threading and len(todo) > 1:
        return _do_calculation_t(
            todo_walk=todo, z=z, force=force, progress=progress)
    else:
        return _do_calculation(
            todo_walk=todo, z=z,
            force=force, progress=progress, verbose=verbose,
        )


def do_all_calculations(
        arange, zrange, nshell, n_component=2, formalism='pn',
        dpath_results=DPATH_RESULTS,
        dpath_sources=DPATH_SOURCES,
        dpath_templates=DPATH_TEMPLATES,
        jmin_even=JMIN_EVEN, jmax_even=JMAX_EVEN,
        jmin_odd=JMIN_ODD, jmax_odd=JMAX_ODD,
        force=False, verbose=False, progress=True, threading=True,
):
    """For each proton number (Z) in zrange and mass number (A) in arange,
    copies the directory structure present in sources directory into the
    results directory and does shell and bat calculations for each interaction.
    :param arange: range of mass numbers (A) for which to do the calculations
    :param zrange: range of proton numbers (Z) for which to do the calculations
    :param nshell: major oscillator shell (0=s, 1=p, 2=sd, ...)
    :param n_component: (1 -> neutrons, 2 -> protons and neutrons)
    :param formalism: ('t' -> isospin, 'pn' -> proton/neutron)
    :param dpath_results: path to the results directory
    :param dpath_sources: path to the sources directory
    :param dpath_templates: path to the templates directory
    :param jmin_even: minimum angular momentum for even A
    :param jmax_even: maximum angular momentum for even A
    :param jmin_odd: minimum angular momentum for odd A
    :param jmax_odd: maximum angular momentum for odd A
    :param force: if true, force redoing calculations even if output files
    already exist
    :param verbose: if true, writes regular output of shell and bat
    calculations to stdout; else output is suppressed and written to files
    :param progress: if true, displays a progress bar showing completion of
    jobs (NOTE: not compatible with verbose option)
    :param threading: if true, calculations are multi-threaded
    """
    zrange = list(filter(lambda z0: z0 >= 1, zrange))
    for z in zrange:
        arange0 = list(filter(lambda a: a >= z, arange))
        make_results_dir(
            a_range=arange0, z=z,
            nshell=nshell, ncomponent=n_component, formalism=formalism,
            jmin_even=jmin_even, jmax_even=jmax_even,
            jmin_odd=jmin_odd, jmax_odd=jmax_odd,
            dpath_results=dpath_results,
            dpath_sources=dpath_sources,
            dpath_templates=dpath_templates,
            force=force,
        )
        make_usdb_dir(
            a_range=arange0, z=z, nshell=nshell,
            jmin_even=jmin_even, jmax_even=jmax_even,
            jmin_odd=jmin_odd, jmax_odd=jmax_odd,
            dpath_results=dpath_results,
            dpath_templates=dpath_templates,
            force=force,
        )
        remove_empty_directories(dpath_results, remove_root=False)
        do_calculations(
            a_range=arange0, z=z, force=force, dpath_results=dpath_results,
            verbose=verbose, progress=progress, threading=threading,
        )


# todo this stuff should be done better
if __name__ == "__main__":
    user_args = argv[1:]
    force0 = False
    verbose0 = False
    while True:
        a0 = user_args[0]
        if a0 == '-f' or a0 == '--force':
            force0 = True
        elif a0 == '-v' or a0 == '--verbose':
            verbose0 = True
        else:
            break
        user_args = user_args[1:]
    if len(user_args) == 5:
        nshell0 = int(user_args[0])
        ncomponent0 = int(user_args[1])
        formalism0 = user_args[2]
        amin, zmin = [int(x) for x in user_args[3:]]
        do_all_calculations(
            arange=[amin], zrange=[zmin],
            nshell=nshell0, n_component=ncomponent0,
            formalism=formalism0, force=force0, verbose=verbose0
        )
    elif len(user_args) == 6:
        nshell0 = int(user_args[0])
        ncomponent0 = int(user_args[1])
        formalism0 = user_args[2]
        amin, amax, zmin = [int(x) for x in user_args[3:]]
        do_all_calculations(
            arange=range(amin, amax+1), zrange=[zmin],
            nshell=nshell0, n_component=ncomponent0,
            formalism=formalism0, force=force0, verbose=verbose0
        )
    elif len(user_args) == 7:
        nshell0 = int(user_args[0])
        ncomponent0 = int(user_args[1])
        formalism0 = user_args[2]
        amin, amax, zmin, zmax = [int(x) for x in user_args[3:]]
        do_all_calculations(
            arange=range(amin, amax+1), zrange=range(zmin, zmax+1),
            nshell=nshell0, n_component=ncomponent0,
            formalism=formalism0, force=force0, verbose=verbose0
        )
    else:
        print ('User entered %d arguments. ' % (len(user_args),) +
               'shell_calc.py requires 5-7 arguments.')
