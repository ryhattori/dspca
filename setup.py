from distutils.core import setup
import setuptools

setup(
    name='dspca',
    version='1.0.1',
    packages=setuptools.find_packages(),
    url='http://github.com/ryhattori/dspca',
    license='MIT License',
    author='Ryoma Hattori',
    author_email='rhattori0204@gmail.com',
    description='Demixed subspace principal component analysis (dsPCA)',
    install_requires=[
        "numpy",
        "scikit-learn",
        "scipy",
        "matplotlib",
    ]
)
