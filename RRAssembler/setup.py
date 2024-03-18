from distutils.core import setup, Extension
from distutils.extension import Extension

from Cython.Build import cythonize
import numpy
import pysam


def make_ext(modname, pyxfilename):
    from distutils.extension import Extension
    import pysam
    return Extension(name=modname,
          sources=[pyxfilename],
          extra_link_args=pysam.get_libraries(),
          include_dirs=pysam.get_include(),
          define_macros=pysam.get_defines())


setup(
    ext_modules=[
        make_ext("recursive_find", "recursive_find.c")
        # Extension("recursive_find", ["recursive_find.c"],
        #           include_dirs=[numpy.get_include(), pysam.get_include()]),
        #           define_macros=pysam.get_defines(),
        #           extra_link_args=pysam.get_libraries()
    ],
)

# Or, if you use cythonize() to make the ext_modules list,
# include_dirs can be passed to setup()

setup(
    ext_modules=cythonize("recursive_find.pyx"),
    include_dirs=[numpy.get_include(), pysam.get_include(), "/usr/local/include/htslib"]
)    

