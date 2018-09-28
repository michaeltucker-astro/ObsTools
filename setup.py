import subprocess, pip, shutil
import os, sys


PACKAGES = [
	'SpectRes',
	'astropy',
	'numpy',
	'matplotlib',
	'argparse'
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

	resp = input('Make Python scripts executables? [Y/n] >').strip().lower()
	if resp == 'y':
		MakeExecutables()

	resp2 = input('Add this directory to $PATH? [Y/n] >').strip().lower()
	if resp == 'y':
		AddToPath()


def CheckPkg(pkg):
	try:
		exec('import %s' % pkg)
	except [ImportError, ModuleNotFoundError]:
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

		subprocess.run('chmod u+x '+ff)


def AddToPath():
	bashrc = '~/.bashrc'
	if not os.path.isfile(bashrc): 
		print('WARNING (AddToPath): Could not find ~/.bashrc file, cannot add to path!')
		return

	cwd = os.getcwd()
	PATHline = 'export PATH="%s:$PATH"\n' % cwd
	for line in open(bashrc,'r').readlines():
		if line == PATHline: return

	with open(bashrc, 'a+') as ofile:
		ofile.write(PATHline)

	subprocess.run('source '+bashrc)
	print('Successfully added to $PATH!')

if __name__=='__main__':
	main()