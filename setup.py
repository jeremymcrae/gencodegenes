
import glob
import os
import sys
import shutil

from setuptools import setup, Extension
from Cython.Build import cythonize

EXTRA_COMPILE_ARGS = []
EXTRA_LINK_ARGS = []
if sys.platform == "win32":
    EXTRA_COMPILE_ARGS += ['/std:c++14']
else:
    EXTRA_COMPILE_ARGS += ['-std=c++11']
    if sys.platform == "darwin":
        EXTRA_COMPILE_ARGS += ["-stdlib=libc++"]
        EXTRA_LINK_ARGS += ["-stdlib=libc++"]
        sdk_base = "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk"
        if os.path.exists(sdk_base):
            EXTRA_COMPILE_ARGS += [
                f"-I{sdk_base}/usr/include/c++/v1",
                f"-I{sdk_base}/usr/include",
            ]
            EXTRA_LINK_ARGS += [
                f"-L{sdk_base}/usr/lib",
            ]

gencode_sources = [
    "src/gencodegenes/gencode.pyx",
    "src/gencode.cpp",
    "src/gtf.cpp",
    "src/tx.cpp",
]

libs = ['z']
include_dirs = ['src/']

if sys.platform == 'win32':
    gencode_sources += glob.glob('src/zlib/*.c')
    include_dirs.append('src/zlib/')
    libs = []

extensions = [
    Extension("gencodegenes.transcript",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=[
            "src/gencodegenes/transcript.pyx",
            "src/tx.cpp"],
        include_dirs=["src/"],
        language="c++"),
    Extension("gencodegenes.gencode",
        extra_compile_args=EXTRA_COMPILE_ARGS,
        extra_link_args=EXTRA_LINK_ARGS,
        sources=gencode_sources,
        include_dirs=include_dirs,
        libraries=libs,
        language="c++"),
    ]

# include tx.h in the package, for downstream usage
shutil.copy("src/tx.h", "src/gencodegenes/tx.h")
shutil.copy("src/tx.cpp", "src/gencodegenes/tx.cpp")

setup(
    package_dir={'': 'src'},
    ext_modules=cythonize(extensions),
    )
