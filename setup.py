import subprocess, pip, shutil
import os, sys

PACKAGES = [
	'SpectRes',
	'astropy',
	'numpy',
	'matplotlib',
	'argparse',
	'pathlib',
	'ipdb',
	'pandas'
]

FileList = [
	'quick-class.py',
	'obstools.py',
	'extractSNIFS.py'
]


def main():
	cwd = os.getcwd()

	for pkg in PACKAGES:
		CheckPkg(pkg)

	CheckDependencies()

	resp = input('Make Python scripts executables? [Y/n] >').strip().lower()
	if resp == 'y' or resp == '':
		MakeExecutables()

	resp2 = input('Add this directory to $PATH? [Y/n] >').strip().lower()
	if resp == 'y' or resp =='':
		AddToPath()


def CheckPkg(pkg):
	try:
		exec('import %s' % pkg)
	except (ImportError, ModuleNotFoundError):
		print('Failed to find package: %s, running pip install' % pkg)
		subprocess.call([sys.executable, "-m", "pip", "install", pkg])
	finally:
		exec('import %s' % pkg)
		return

def MakeExecutables():
	python_cmd = shutil.which('python')
	startline = '#!'+python_cmd+'\n'
	for ff in FileList:
		lines = open(ff, 'r').readlines()
		if lines[0].startswith('#!'): continue
		lines = [startline]+lines
		with open(ff, 'w') as ofile:
			ofile.write(''.join(lines))


		subprocess.run('chmod +x '+ff, shell=True)


def AddToPath():
	import pathlib
	import ipdb

	bashrc = str(pathlib.Path.home())+'/.bashrc'
	if not os.path.exists(bashrc): 
		print('WARNING (AddToPath): Could not find %s file, cannot add to path!' % bashrc)
		return

	cwd = os.getcwd()
	PATHline = 'export PATH="%s:$PATH"\n' % cwd
	for line in open(bashrc,'r').readlines():
		if line == PATHline: return

	with open(bashrc, 'a+') as ofile:
		ofile.write(PATHline)

	print('Successfully added to $PATH! ')
	print('Run "source %s" to update current shell PATH' % bashrc)

def CheckDependencies():
	dependents = [
	'google-chrome'
	]

	for dep in dependents:
		path = shutil.which(dep)
		if path == None:
			print('WARNING (CheckDependencies): could not detect %s executable, edit FinderChart (obstools.py) accordingly!' % dep)

if __name__=='__main__':
	main()