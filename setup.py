from setuptools import setup
from setuptools_rust import Binding, RustExtension

setup(
    name="Gemi3",
    version="1.0",
    rust_extensions=[RustExtension("Gemi3.rsCAVI", binding=Binding.PyO3)],
    packages=["Gemi3","Gemi3.pyCAVI"],
    # rust extensions are not zip safe, just like C-extensions.
    zip_safe=False,
)