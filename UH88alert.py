from datetime import datetime
import numpy as np
from twilio.rest import Client as SMSclient
import requests, os, time
import argparse

#Twilio account tokens, in account SID/ auth-token pairs
TOKENS = {
	'MAT':{'auth-token':'daa9733bef28399e1efaa38f61920ebd',
			'accountSID':'AC6eaf6e92547d83a80ad9f3920f2dc324',
			'actual_phone':'+19197104630',
			'twilio_phone':'+12526741988',
			'email':'tuckerma@hawaii.edu'
		}
}

WeatherURL = 'http://uh88.ifa.hawaii.edu/weather/'

CRITICAL_HUMIDITY = 85.0 #percent
WARN_HUMIDITY = 70.0 #percent

CRITICAL_WIND_SPEED = 50.0 #mph
WARN_WIND_SPEED = 35.0 #mph

CRITICAL_TEMPERATURE = 20.0 #Celsius
WARN_TEMPERATURE = 10.0 #Celsius

def main(initials, atype, update=0, interval=5.0, start=0.0, end=6.0, logfile='UH88alerts.log', verbose=False):

	if verbose:
		print('#'*30)
		print('UH88 weather alert system')
		print('Selected setup:')
		print('\tInitials: %s' % initials)
		print('\tAlert type: %s' % atype)
		print('\tWait to start: %3g hr' % start)
		print('\tTracking duration: %3g hr' % end)
		print('\tWeather refresh interval: %3g minutes' % interval)
		print('\tUpdate interval: %3g minutes' % update)
		print('\tLog file: %s\n' % logfile)

	if initials not in TOKENS.keys(): raise ValueError('Initials %s not in TOKEN database!' % initials)
	if start >= 0.0:
		if verbose: print('Executing sleep for %3g hr at %s' % (start, datetime.now()))
		time.sleep(3600.0*start)

	if verbose: print('Starting UH88 alert system')
	alerts = Alert(initials, atype, interval, update, logfile, end, verbose)
	alerts.LOOP()




class Alert():
	def __init__(self, initials, atype='text', interval=5.0, update=0, logfile='UH88alerts.log', duration=6.0, verbose=False):
		self.initials = initials
		if atype == 'both':
			self.email = True
			self.text = True
		elif atype == 'email':
			self.email=True
			self.text = False
		elif atype == 'text':
			self.email=False
			self.text = True
		elif atype == 'none':
			self.email=False
			self.text = False
		else:
			raise RuntimeError('Unrecognized alert type: %s' % atype)

		if verbose: print('Opening log file.')
		self.logger = Logger(logfile)


		if update != 0 and update < interval: 
			if verbose: print('WARNING: update < interval, setting update = interval')
			self.logger('update: %lf min, interval: %lf min' % (update, interval))
			self.logger('update < interval, setting update = interval')

		self.interval = interval
		self.update = update
		self.duration = duration


		self.verbose = verbose

		self.humidity_perc = []
		self.temperature = []
		self.windspeed = []

		self.running = True
		self.WARN = False
		self.ALARM = False

		self.time_start = datetime.now()
		self.last_update = datetime.now()

		self.logger('Starting UH88 alert system: %s' % str(datetime.now()))
		if self.verbose: print('Alert system primed and ready')

	def LOOP(self):
		if self.verbose: print('\tStarting main loop')
		while self.running:
			self.Update()
			self.log()
			self.CheckTime()
			self.wait()

	def Update(self):
		if self.verbose: print('Running Update at %s' % str(datetime.now()))

		self.humidity, self.temp, self.wind_speed = ScrapeWeatherPage(5)

		if self.verbose: print('\tCurrent: Humidity: %3.1f%% Temp: %3.1f C WindSpeed: %3.1f mph' % (self.humidity, self.temp, self.wind_speed))

		if ( (self.humidity >= CRITICAL_HUMIDITY) or (self.temp >= CRITICAL_TEMPERATURE) or (self.wind_speed >= CRITICAL_WIND_SPEED) ):
			if self.ALARM: pass
			else:  
				self.SendAlert('alarm')
				self.ALARM = True

		elif ( (self.humidity >= WARN_HUMIDITY) or (self.temp >= WARN_TEMPERATURE) or (self.wind_speed >= WARN_WIND_SPEED) ):
			if self.ALARM:
				self.SendAlert('warn')
				self.ALARM = False

			elif not self.WARN and not self.ALARM:
				self.SendAlert('warn')
				self.WARN = True

		else:
			if self.ALARM or self.WARN:
				self.SendAlert('safe')
				self.ALARM = False
				self.WARN = False

		self.humidity_perc.append(self.humidity)
		self.temperature.append(self.temp)
		self.windspeed.append(self.wind_speed)

	def CheckTime(self):
		now = datetime.now()
		elapsed = (now - self.time_start).total_seconds()/3600.0
		if self.verbose: print('Current elapsed time: %3g hr' % elapsed)
		if elapsed > self.duration:self.Close()

		if self.update == 0.0: return
		elapsed = (now - self.last_update).total_seconds()/60.0
		if self.verbose: print('\nElapsed time since last update: %3g minutes' % elapsed)
		if elapsed >= self.update: 
			if self.verbose: print('\tCriteria met, Sending update!')
			self.SendUpdate()
			self.last_update = datetime.now()

	def log(self):
		self.logger('='*30)
		self.logger('DATETIME: %s' % str(datetime.now()))
		self.logger('WARN: %r' % self.WARN)
		self.logger('ALARM: %r' % self.ALARM)
		self.logger('Humidity: %3.1f%%' % self.humidity)
		self.logger('Temperature: %3.1f C' % self.temp)
		self.logger('Wind speed: %3.1f mph\n' % self.wind_speed)

	def SendAlert(self, message):
		if self.verbose: print('Sending alert!')
		if message == 'safe':
			body = 'UH88 Weather Update: Current Status SAFE'
		elif message == 'warn':
			body = 'UH88 Weather Update: Current status WARN'
		elif message == 'alarm':
			body = 'UH88 Weather Update: Current Status CRITICAL'
		else:
			raise ValueError('Unknown message: %s' % message)

		body += '\nHumidity: %3.1f%%\nTemp: %3.1fC\nWind: %3.1fmph' % (self.humidity, self.temp, self.wind_speed)

		if self.text:
			SendTextAlert(self.initials, body)
		if self.email:
			SendEmailAlert(self.initials, body)

	def SendUpdate(self):
		body = 'UH88 Weather Update:\nHumidity: %3.1f%%\nTemp: %3.1fC\nWind speed: %3.1fmph\n' % (self.humidity, self.temp, self.wind_speed)
		if self.email:
			SendEmailAlert(self.initials, body)
		if self.text:
			SendTextAlert(self.initials, body)

	def wait(self):
		time.sleep(60.0*self.interval)

	def Close(self):
		if self.verbose: print('\nClose called, exiting alert system...\n')
		self.logger('Exiting alert system: %s' % str(datetime.now()))
		self.logger.close()
		self.running = False

class Logger():
	def __init__(self, fname):
		self.fname = fname
		resp = 'null'
		while os.path.exists(self.fname) or resp.strip().lower() not in ['', 'y', 'n', 'yes', 'no']:
			resp = input('Log file %s already exists, overwrite? [Y/n] >' % self.fname)
			if resp.strip().lower() in ['n', 'no']: 
				self.fname = input('\tNew filename >').strip()
			elif resp.strip().lower() in ['', 'y', 'yes']:
				break
			else:
				print('\tUnrecognized response: %s' % resp.strip().lower())

		self.file = open(self.fname, 'w')

	def __call__(self, string):
		self.file.write(string+'\n')

	def close(self):
		self.file.close()

def SendTextAlert(initials, body):

	info = TOKENS[initials] 
	client = SMSclient(info['accountSID'], info['auth-token'])
	sender, recipient = [TOKENS[initials][key] for key in ['twilio_phone', 'actual_phone']]
	message = client.messages.create(
		from_=sender,
		body=body,
		to=recipient
		) 

def SendEmailAlert(initials, body):
	email_addr = TOKENS[initials]['email']


def ScrapeWeatherPage(retries=1, wait=60.0):
	resp = requests.get(WeatherURL)
	if resp.status_code != 200:
		if retries == 0: return None
		print('WARNING: Bad response from Weather page, conducting %d retries...' % retries)
		for attempt in range(retries):
			time.wait(wait)
			results = ScrapeWeatherPage(0)
			if type(results) == tuple: return results
		raise IOError('Failed to retrieve data from UH88 weather page, check internet connection!')

	data = resp.content.decode().splitlines()[31].strip().split()
	temp = float(data[2].strip()[:-1])
	humid = float(data[5].strip()[:-1])
	wind = float(data[-1].strip()[:-4])
	return humid, temp, wind

if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Sends text/email alerts about Uh88 weather conditions.')
	parser.add_argument('--interval', '-i', help='Interval between weather data refresh (minutes). Default: 5.0', default=5.0, type=float)
	parser.add_argument('--person', '-p', help='Initials of who to send alerts to. Default: MAT', default='MAT', type=str, 
						choices = list(TOKENS.keys()))
	parser.add_argument('--atype', '-a', help='Alert type: email, text, or both. Default: text', default='text', type=str, 
						choices=['email', 'text', 'both'])
	parser.add_argument('--log', '-l', help='Log file. Default: UH88alerts.log', default='UH88alerts.log', type=str)
	parser.add_argument('--update', '-u', help='Interval between weather updates sent [min] (0 = off). Default: 0', default=0.0, type=float)
	parser.add_argument('--start', '-s', help='Delay between now and when to start tracking weather (hours). Default: 0', default=0.0, type=float)
	parser.add_argument('--end', '-e', help='How long to track weather (hours). Default: 6.0', default=1.0, type=float)
	parser.add_argument('--verbose', '-v', help='Verbose? [0/1] Default: 0', default=0, type=int, choices=[0,1])
	args = parser.parse_args()

	main(args.person, args.atype, args.update, args.interval, args.start, args.end, args.log, args.verbose)
