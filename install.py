#!/usr/bin/env python
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# hap.py installer:
#
# * checks dependencies
# * builds code
# * makes virtualenv

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
import glob
import fnmatch
import multiprocessing
import urllib2


def check_python_version():
    """ Check if the python version is sufficient """
    if sys.version_info < (2, 7, 3):
        raise Exception("You will need to run this with Python >= 2.7.3")


def create_python_environment(source_dir, args):
    """ Create a Python runtime environment
    :return: shebang with path to the python executable
    """
    interp = args.python_interp
    pyver = subprocess.check_output(interp + " -c \"import sys; print ','.join(map(str, list(sys.version_info["
                                    "0:3])))\"", shell=True).strip().split(",")
    pyver = tuple(map(int, pyver))

    if pyver < (2, 7, 3):
        raise Exception("Python >= 2.7.3 is required for installation.")

    # system python -- just return interp
    if args.python == "system":
        return "#!" + interp

    if not args.python_venv_dir:
        raise Exception("Please specify a virtualenv target installation directory.")

    if args.python_venv_dir_force:
        try:
            shutil.rmtree(args.python_venv_dir)
        except:
            pass

    if os.path.exists(args.python_venv_dir) and not args.python_venv_dir_update:
        raise Exception("The virtual environment directory already exists.")

    virtualenv_tempdir = tempfile.mkdtemp(prefix="virtualenv", dir=args.scratch_path)
    try:
        ve_tgz = os.path.join(source_dir, "external", "virtualenv-12.0.7.tar.gz")
        to_run = "cd %s && tar xzf %s" % (virtualenv_tempdir, ve_tgz)
        print >>sys.stderr, to_run
        subprocess.check_call(to_run, shell=True)

        ve_exec = os.path.join(virtualenv_tempdir, "virtualenv-12.0.7", "virtualenv.py")
        to_run = "%s -p %s %s" % (ve_exec, interp, args.python_venv_dir)
        print >>sys.stderr, to_run
        subprocess.check_call(to_run, shell=True)
    finally:
        if not args.keep_scratch:
            try:
                shutil.rmtree(virtualenv_tempdir)
            except:
                pass

    # install requirements
    ve_python = os.path.join(args.python_venv_dir, "bin", "python")
    ve_pip = os.path.join(args.python_venv_dir, "bin", "pip")

    deleteme = None
    try:
        cmds = [ve_pip, "install", "--no-cache-dir"]

        if args.fix_cert:
            response = urllib2.urlopen('http://curl.haxx.se/ca/cacert.pem')
            certdata = response.read()
            f = tempfile.NamedTemporaryFile(delete=False)
            deleteme = f.name
            f.write(certdata)
            f.close()
            cmds.insert(1, " --cert")
            cmds.insert(2, deleteme)

        for x in open(os.path.join(source_dir, "happy.requirements.txt")):
            print >>sys.stderr, " ".join(cmds + [x])
            subprocess.check_call(" ".join(cmds + [x]), shell=True)
    finally:
        if deleteme:
            os.unlink(deleteme)

    return "#!" + ve_python


def replace_shebang(filename, shebang):
    """ Replace shebang line / reheader script files """
    print >>sys.stderr, "Fixing shebang line in " + filename

    with open(filename) as f:
        lines = f.readlines()

    with open(filename, "w") as f:
        removed = False
        print >> f, shebang
        for i, l in enumerate(lines):
            if not removed and l.startswith("#!") and i < 10:
                removed = True
            else:
                f.write(l)


def build_haplotypes(source_dir, build_dir, args):
    if args.boost:
        boost_prefix = "BOOST_ROOT=%s " % args.boost
    else:
        boost_prefix = ""
    config_command = "%s/configure.sh %s %s %s" % (source_dir, args.configuration, args.setup,
                                                   args.targetdir)
    if args.sge:
        config_command += " -DUSE_SGE=ON"

    if args.build_rtgtools:
        config_command += " -DBUILD_VCFEVAL=ON"
        if args.rtgtools_wrapper:
            if not os.path.exists(args.rtgtools_wrapper):
                raise Exception("RTG-tools wrapper %s doesn't exist." % args.rtgtools_wrapper)
            config_command += "-DVCFEVAL_WRAPPER=%s" % \
                              os.path.abspath(args.rtgtools_wrapper).replace(" ", "\\ ")

    to_run = boost_prefix + "cd %s && %s %s" % (build_dir, boost_prefix, config_command)
    print >>sys.stderr, to_run
    subprocess.check_call(to_run, shell=True)

    setupscript = ""
    if args.setup != "auto":
        setupscript = " . %s && " % os.path.join(source_dir, "src", "sh", args.setup + "-setup.sh")

    setupscript += boost_prefix

    to_run = setupscript + "cd %s && %s make -j%i" % (build_dir, setupscript, args.processes)
    print >>sys.stderr, to_run
    subprocess.check_call(to_run, shell=True)

    to_run = setupscript + "cd %s && %s make -j%i install" % (build_dir, setupscript, args.processes)
    print >>sys.stderr, to_run
    subprocess.check_call(to_run, shell=True)


def test_haplotypes(source_dir, python_shebang, args):
    """ Run the unit + integration tests
    """
    to_run = "cd %s && %s" % (args.targetdir, os.path.join(source_dir, "src", "sh", "run_tests.sh"))
    print >>sys.stderr, to_run
    os.environ["PYTHON"] = python_shebang[2:]
    subprocess.check_call(to_run, shell=True)


def main():
    check_python_version()

    source_dir = os.path.abspath(os.path.dirname(__file__))

    parser = argparse.ArgumentParser("hap.py installer")
    parser.add_argument("targetdir", help="Target installation directory")

    parser.add_argument("--sge-mode",
                        dest="sge",
                        action="store_true",
                        default=False,
                        help="Enable SGE mode, which will require an additional command "
                             "line option \"--force-interactive\" to run interactively.")

    parser.add_argument("--python",
                        dest="python",
                        choices=["system", "virtualenv"],
                        default="system",
                        help="Which Python to use in the installation. 'virtualenv' "
                             "will create a virtual environment in the folder "
                             "specified with --python-virtualenv-dir")

    parser.add_argument("--python-interpreter", dest="python_interp",
                        default=sys.executable,
                        help="Python interpreter to use for the installed hap.py.")

    parser.add_argument("--python-virtualenv-update", dest="python_venv_dir_update",
                        default=False, action="store_true",
                        help="Update virtualenv if it already exists.")

    parser.add_argument("--python-virtualenv-force", dest="python_venv_dir_force",
                        default=False, action="store_true",
                        help="Force creating a virtualenv even if the target directory"
                             " already exists. USE WITH CARE, THIS WILL REMOVE THE "
                             "VIRTUALENV DIRECTORY!")

    parser.add_argument("--python-virtualenv-dir", dest="python_venv_dir", default="",
                        help="Directory to install the virtualenv in.")

    parser.add_argument("--pip-fix-cert", dest="fix_cert", default=False, action="store_true",
                        help="Download and use certificate file in case of Linux distributions"
                             " which have an outdated certificate file which makes pip fail.")

    # C++ compile options
    setups = [os.path.basename(x).replace("-setup.sh", "")
              for x in glob.glob(os.path.join(source_dir, "src", "sh", "*-setup.sh"))]

    setups.insert(0, "auto")

    parser.add_argument("--configuration", dest="configuration",
                        choices=["Debug", "Release", "RelWithDebInfo", "install"],
                        default="Release", help="Build configuration (use Release if unsure).")

    parser.add_argument("--setup", dest="setup", choices=setups,
                        default="auto", help="Build setup (or auto to use system-wide packages).")

    parser.add_argument("--boost-root", dest="boost",
                        help="Where to find Boost.", default="")

    parser.add_argument("--scratch-path", dest="scratch_path",
                        help="Where to build.", default="/tmp")

    parser.add_argument("--keep-scratch", dest="keep_scratch",
                        help="Keep the scratch folder.", default=False,
                        action="store_true")

    parser.add_argument("--with-rtgtools", dest="build_rtgtools", default=False,
                        action="store_true",
                        help="Get and build rtgtools. You need to have Java and Ant for this.")

    parser.add_argument("--rtgtools-wrapper", dest="rtgtools_wrapper", default=None,
                        help="Wrapper script for rtgtools. This is optional, it is useful "
                             "when the default version of Java must be replaced / the environment "
                             "needs updating. There is an example in src/sh/rtg-wrapper.sh.")

    parser.add_argument("--build-processes", dest="processes",
                        default=multiprocessing.cpu_count(), type=int,
                        help="Number of parallel processes to use for building.")

    parser.add_argument("--no-tests", dest="run_tests", default=True, action="store_false",
                        help="Disable unit tests")

    args = parser.parse_args()

    args.targetdir = os.path.abspath(args.targetdir)

    if args.python == "virtualenv" and not args.python_venv_dir:
        raise Exception("Please specify a virtualenv target installation directory.")

    if args.python_venv_dir and not args.python_venv_dir.startswith("/"):
        args.python_venv_dir = args.targetdir

    if "LD_LIBRARY_PATH" in os.environ or "DYLD_LIBRARY_PATH" in os.environ:
        print >>sys.stderr, "WARNING: You have (DY)LD_LIBRARY_PATH set. Make sure these libraries are accessible " \
                            "in the same environment you will run in."

    # fix dynamic linking
    if "LD_LIBRARY_PATH" in os.environ:
        os.environ["LD_RUN_PATH"] = os.environ["LD_LIBRARY_PATH"]

    # check/make Python environment
    python_shebang = create_python_environment(source_dir, args)

    if args.boost and not os.path.exists(args.boost):
        raise Exception("Boost directory doesn't exist.")

    # build hap.py
    build_dir = tempfile.mkdtemp(prefix="build", dir=args.scratch_path)
    try:
        build_haplotypes(source_dir, build_dir, args)
    finally:
        if not args.keep_scratch:
            shutil.rmtree(build_dir)

    # reheader Python files
    for root, _, filenames in os.walk(args.targetdir):
        for filename in fnmatch.filter(filenames, '*.py'):
            replace_shebang(os.path.join(root, filename), python_shebang)

    if args.run_tests:
        test_haplotypes(source_dir, python_shebang, args)


if __name__ == '__main__':
    main()
