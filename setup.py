from distutils.core import setup

# Dynamically calculate the version based on pymbar.VERSION.
version = __import__('pymbar').__version__

setup(
    name = "pymbar",
    version = version.replace(' ', '-'),
    url = 'https://simtk.org/home/pymbar',
    download_url = 'http://github.com/davecap/pymbar',
    author = 'Shirts MR and Chodera JD',
    author_email = 'mrshirts@gmail.com',
    description = 'A Python implementation of the multistate Bennett acceptance ratio (MBAR)',
    packages = ['pymbar'],
)

