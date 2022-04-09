import glob
import os
import platform
from waflib.TaskGen import extension, feature, task_gen
from waflib.Task import Task
from waflib import Build

VERSION = '0.0'
APPNAME = 'bvh8tool'

top = '.'


def options(opt):
    if ('COMSPEC' in os.environ):
        opt.load('msvc')
    else:
        opt.load('clang++')

def configure(conf):
    if ('COMSPEC' in os.environ):
        conf.load('msvc')
    else:
        conf.load('clang++')

def build(bld):

    if ('COMSPEC' in os.environ):
        platform_defines = ['_CRT_SECURE_NO_WARNINGS', 'PLATFORM_WINDOWS']
    else:
        platform_defines = ['_CRT_SECURE_NO_WARNINGS', 'PLATFORM_LINUX']
    includes = ['source', 'includes']

    # RELEASE
    sdk_lib_path = []
    libs = []
    compile_flags = ['-O0', '-g', '-msse4.1']
    linker_flags = []

    # Build risctool
    bld.program(
        source=glob.glob('*.cpp'),
        cxxflags=compile_flags,
        ldflags=linker_flags,
        target='bvh8tool',
        defines=platform_defines,
        includes=includes,
        libpath=[sdk_lib_path],
        lib=libs)
