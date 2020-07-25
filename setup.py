from setuptools import setup

setup(
    name="KerrNewman",
    version="0.1",
    author='Gabriele Bozzola',
    author_email='gabrielebozzola@arizona.edu',
    packages=['KerrNewman'],
    scripts=['bin/estimate_spins.py'],
    description='Compute ISCO of Kerr-Newman black holes and estimate final spin of mergers',
    install_requires=["numpy", "scipy", "sympy"],
    license='LICENSE.txt',
)
