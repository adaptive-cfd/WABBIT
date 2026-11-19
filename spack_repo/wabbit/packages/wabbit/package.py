# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

# ----------------------------------------------------------------------------
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install wabbit
#
# You can edit this file again by typing:
#
#     spack edit wabbit
#
# See the Spack documentation for more information on packaging.
# ----------------------------------------------------------------------------

from spack_repo.builtin.build_systems.makefile import MakefilePackage
from spack.package import *


class Wabbit(MakefilePackage):
    """Wavelet Adaptive Block-Based solver for Interactions with Turbulence"""

    homepage = "www.cfd.tu-berlin.de/"
    url = "https://github.com/adaptive-cfd/WABBIT/archive/refs/tags/v1.0-CiCP21.tar.gz"

    maintainers("tommy-engels", "Philipp137")

    license("GPL-3.0", checked_by="freifrauvonbleifrei")

    version("master", branch="master")
    version(
        "1.0-CiCP21", sha256="212ac891f4d97a42a0c43b38e2bc867f18628045c6353f64774190fce6e75a88"
    )

    # variant("test", default=False, description="Build WABBIT with tests")
    # #TODO git clone https://github.com/adaptive-cfd/python-tools.git

    depends_on("fortran", type="build")

    depends_on("blas", type=("build", "link"))
    depends_on("fftw", type=("link"), when="@2.0:")
    depends_on("hdf5+mpi+fortran", type=("build", "run"))
    depends_on("lapack", type=("build", "link"))
    depends_on("mpi")
    depends_on("py-h5py", type=("run"))
    depends_on("py-numpy", type=("run"))

    def setup_build_environment(self, env: EnvironmentModifications) -> None:
        env.set("WABBIT_BLAS_ROOT", self.spec["blas"].prefix.lib)
        env.set("WABBIT_LAPACK_ROOT", self.spec["lapack"].prefix.lib)
        if "^fftw" in self.spec:
            env.set("FFT_ROOT", self.spec["fftw"].prefix)
        env.set("HDF_ROOT", self.spec["hdf5"].prefix)
        if "^hdf5@1.13:" in self.spec:
            env.set("HDF_SOURCE", self.spec["hdf5"].prefix.include)
        # if test, set WABBIT_PT_ROOT=...

    def install(self, spec, prefix):
        mkdir(prefix.bin)
        install("wabbit", prefix.bin)
        # install_tree("libwabbit.a", prefix.lib)

    def edit(self, spec, prefix):
        makefile = FileFilter("LIB/fortran.mk")
        makefile.filter("blas", self.spec["blas"].name)
        makefile.filter("lapack", self.spec["lapack"].name)
