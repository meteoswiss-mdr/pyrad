from setuptools import setup

setup(
    name='DualPol',
    version='1.0.2',
    author='Timothy Lang',
    author_email='timothy.j.lang@nasa.gov',
    packages=['dualpol', ],
    license='LICENSE.md',
    description='Python Interface to Dual-Pol Radar Algorithms (DualPol)',
    long_description=open('description.txt').read(),
    install_requires=['arm_pyart', 'csu_radartools', 'skewt'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console"
        ],
)
