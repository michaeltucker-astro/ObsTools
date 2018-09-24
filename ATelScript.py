import os
from astropy.time import Time

def WriteATel(outfile='ATel.txt'):

	Names = []
	IAUnames = []
	redshift = []
	zsrc = []
	OBSdate = []
	dtype = []
	phases = []
	DiscATel = []
	redshift_src = []
	Notes = []
	NoteValues = {}

	while True:
		name = input('\nName of object,IAU designation (q to quit)>')
		if name == 'q': break
		Names.append(name.split(',')[0])
		IAUnames.append(name.split(',')[1])
		z = input('Redshift (src) > ').split(',')
		redshift.append(z)
		date = Time(float(input('MJD obs >')), format='mjd')
		date = str(date.datetime.year)+'-'+str(date.datetime.month)+'-'+str(date.datetime.day)+(str(round(date.datetime.hour/24.0,2)))[1:]
		OBSdate.append(date)
		dtyp, phase = input('Type,Phase > ').split(',')
		dtype.append(dtype)
		phases.append(phase)
		atel = input('Disc. ATel #')
		if atel != '':
			DiscATel.append(atel)
		else: DiscATel.append('')

		note = input('Notes (csv) >')
		if note == '': Notes.append('')
		else: Notes.append(note)

	Nobj = len(Names)
	print('Received %d classified objects.' % Nobj)

	header = [
	'The Spectral Classification of Astronomical Transients (SCAT) survey (ATel #11444) presents the classification of %d optical transients.' % Nobj,
	'We report optical spectroscopy (330-970nm) taken with the University of Hawaii 88-inch (UH88) telescope using the SuperNova Integral Field Spectrograph (SNIFS).',
	'Transients were classified using the SuperNova IDentification code (SNID, Blondin & Tonry 2007, ApJ, 666, 1024).',
	'\n\nSurvey Name | IAU Name | Date Obs. | Disc ATel. | Type | Phase | Redshift | Notes\n',
	'-'*90
	]

	ofile = open(outfile, 'w')
	ofile.write(' '.join(header))
	print(Names, IAUnames, OBSdate, redshift, dtype, phases, DiscATel, Notes)
	for item in zip(Names, IAUnames, OBSdate, DiscATel, dtype, phases, redshift, Notes):
		line = '%s\t%s\t%s\t%s\tATel #%s\t%s (%s)\t%s\t%s\t%s\n' % item
		print (line)
		ofile.write(line)
	ofile.close()