
import os
import subprocess
import sys
import shutil

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
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

ZLIB_NG_DIR = os.path.abspath('src/zlib-ng')

class build_ext_zlib_ng(build_ext):
    ''' build zlib-ng as a static library (via cmake), before building extensions
    
    zlib-ng is built with its native API (zng_ prefixed functions), so it can't
    clash with any system zlib loaded into the python process.
    '''
    def run(self):
        build_dir = os.path.abspath(os.path.join(self.build_temp, 'zlib-ng'))
        config = [
            '-DCMAKE_BUILD_TYPE=Release',
            '-DCMAKE_POSITION_INDEPENDENT_CODE=ON',
            '-DBUILD_SHARED_LIBS=OFF',
            '-DZLIB_COMPAT=OFF',
            '-DWITH_GZFILEOP=ON',
            '-DBUILD_TESTING=OFF',
        ]
        archflags = os.environ.get('ARCHFLAGS', '')
        if sys.platform == 'darwin' and archflags:
            # match the architectures requested for the extension (e.g. by cibuildwheel)
            archs = ';'.join(archflags.split()[1::2])
            config.append(f'-DCMAKE_OSX_ARCHITECTURES={archs}')
        subprocess.check_call(['cmake', '-S', ZLIB_NG_DIR, '-B', build_dir] + config)
        subprocess.check_call(['cmake', '--build', build_dir, '--config', 'Release',
            '--parallel'])
        
        lib_dir = build_dir
        lib = 'z-ng'
        if sys.platform == 'win32':
            lib_dir = os.path.join(build_dir, 'Release')
            lib = 'zlibstatic-ng'
        for ext in self.extensions:
            if ext.name == 'gencodegenes.gencode':
                ext.include_dirs.append(build_dir)
                ext.library_dirs.append(lib_dir)
                ext.libraries.append(lib)
        super().run()

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
        sources=[
            "src/gencodegenes/gencode.pyx",
            "src/gencode.cpp",
            "src/gtf.cpp",
            "src/tx.cpp"],
        include_dirs=["src/"],
        define_macros=[('WITH_GZFILEOP', None)],
        language="c++"),
    ]

# include tx.h in the package, for downstream usage
shutil.copy("src/tx.h", "src/gencodegenes/tx.h")
shutil.copy("src/tx.cpp", "src/gencodegenes/tx.cpp")

setup(
    package_dir={'': 'src'},
    ext_modules=cythonize(extensions),
    cmdclass={'build_ext': build_ext_zlib_ng},
    )
