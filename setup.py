from distutils.core import setup
import setuptools

with open("README.md", "r") as fh:
    readme = fh.read()

setup(
    name='dspca',
    version='1.0.3',
    packages=setuptools.find_packages(),
    url='http://github.com/ryhattori/dspca',
    license='MIT License',
    author='Ryoma Hattori',
    author_email='rhattori0204@gmail.com',
    description='Demixed subspace principal component analysis (dsPCA)',
    long_description=readme,
    long_description_content_type="text/markdown",
    install_requires=[
        "numpy",
        "scikit-learn",
        "scipy",
        "matplotlib",
    ]
)
