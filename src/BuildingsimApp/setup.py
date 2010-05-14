from distutils.core import setup
import py2exe
import glob

opts = {'py2exe': { 'includes': 'matplotlib.numerix.random_array',
                    'excludes': ['_gtkagg', '_tkagg'],
                    'dll_excludes': ['libgdk-win32-2.0-0.dll','libgobject-2.0-0.dll']
                  }
        }

setup(
    data_files = [(r'matplotlibdata', glob.glob(r'c:\python24\share\matplotlib\*')),
                  (r'matplotlibdata', [r'c:\python24\share\matplotlib\.matplotlibrc'])],
    name = 'buildingsimApp',
    description = 'Simulacion sobre Edificios',
    console=['buildingsimApp.py'])