from disutils.spawn import find_executable
import subprocess, pip
import os, sys


PACKAGES = [
	'SpectRes',
	'astropy',
	'numpy',
	'matplotlib',
	'argparse',
]


def main():
	cwd = os.getcwd()

	for pkg in PACKAGES:
		CheckPkg(pkg)

	resp = input('Make Python scripts executables? [Y/n] >')


def CheckPkg(pkg):
	try:
		exec('import %s' % pkg)
	except [ImportError, ModuleNotFoundError]:
		print('Failed to find package: %s; running pip install %s' % (pkg, pkg))
		pip.main('install', pkg)

	finally:
		exec('import %s' % pkg)
		return
