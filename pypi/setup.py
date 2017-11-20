from setuptools import setup

setup(name='LIMBR',
      version='0.1',
      description='Learning and Imputation for Mass-spec Bias Reduction',
      url='https://github.com/aleccrowell/LIMBR',
      download_url='https://github.com/aleccrowell/LIMBR/releases/tag/0.1'
      author='Alec Crowell',
      author_email='alexander.m.crowell@gmail.com',
      license='MIT',
      keywords=['SVA','SVD','mass-spec','bioinformatics],
      install_requires=['numpy','pandas','scipy','sklearn','itertools','statsmodels','tqdm','pickle','multiprocess','functools'],
      zip_safe=False)
