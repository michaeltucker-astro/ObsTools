import numpy as np
from astropy.io import fits
import requests, json


API_KEY = "vrwmqnznsdwljcsy"


def main():

	login = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": API_KEY})})
	sessionID = login['session']
	
	
