from setuptools import setup, find_packages


def get_requires():
    requires = []
    with open("requirements.txt", "rt") as req_file:
        for line in req_file:
            requires.append(line.rstrip())
    return requires


setup(name='opentargets_pharmgkb',
      version='0.0.1',
      packages=find_packages(),
      install_requires=get_requires(),
      tests_require=get_requires(),
      setup_requires=get_requires(),
      test_suite='tests'
      )
