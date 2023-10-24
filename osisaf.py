import numpy as np
from datetime import date, datetime, timedelta
import glob
import csv
import sys
import os
import urllib3
import shutil
import urllib.request
from contextlib import closing
from math import sqrt, sin, cos, pi, floor, isnan, nan, atan, atan2, ceil
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from PIL import Image, ImageGrab, ImageDraw, ImageFont
import matplotlib.ticker as ticker
import dropbox
from time import sleep

from flask import Flask, request, jsonify, send_file 

monthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
imageFolder = 'https://osisaf.met.no/quicklooks/prod/ice/'

app = Flask(__name__)
@app.route("/")
def hello():
	return "Generate an average sea ice drift map by running a query like this: https://osisaf.onrender.com/average?start=2023-10-12&end=2023-10-16"
@app.get("/average")
def get_average():
#@app.route("/")
#def average():
	print('inside avg')
	start=request.args.get('start')
	end=request.args.get('end')
	if start == None:
		return "Please insert start date (start=YYYY-MM-DD)", 400
	if end == None:
		return "Please insert end date (end=YYYY-MM-DD)", 400
	print(start)
	print(end)
	try:
		parsedstart = datetime.strptime(start,'%Y-%m-%d')
		parsedend = datetime.strptime(end,'%Y-%m-%d')
	except ValueError:
		return "Invalid start of end date. They need to have the format YYYY-MM-DD"
	#if parsedstart == None:
	#	return "Start date should have the format yyyy-mm-dd", 400
	#if parsedend == None:
	#	return "End date should have the format yyyy-mm-dd", 400
	if(parsedstart < datetime(2019,6,1)):
		return "Start date cannot be earlier than 2019-06-01", 400
	if(parsedend-parsedstart).days < 2:
		return "End date should be at least 2 days later than start date", 400
	yesterday = datetime.today() - timedelta(days = 1)
	yesterday = datetime(yesterday.year, yesterday.month, yesterday.day)		
		
	if(parsedend > yesterday):
		return "End date cannot be later than " + dateString(yesterday, "-"), 400
	found = False
	counter = 0
	while(not found and parsedend > yesterday - timedelta(days = 3)):
		try:
			dxx,dyy = loadSimpleFiles(parsedend, parsedend)
			found = True			
		except OSError as e:
			counter += 1
			parsedend = parsedend - timedelta(days = 1);
	if(parsedend-parsedstart).days < 2:
		return "Latest available date is " + dateString(parsedend, "-"), 400
		
	prepareImage(downloadImage(parsedend))
	#yesterday = datetime.today() - timedelta(days = 1)
	#enddate = datetime(yesterday.year, yesterday.month, yesterday.day)
	try:
	     filename = fastSum(parsedstart + timedelta(days = 2),parsedend)
	except OSError as e:
		print(type(e))    # the exception type
		print(e.args)     # arguments stored in .args
		print(e) 
		return("error: " + str(e).replace("data/fast/", "").replace("dx_", ""))
	print('The value of __name__ is ' + __name__)
	return send_file(filename, mimetype='image/png')
	
def downloadImage(date):
	previousDate = date - timedelta(days = 2)
	filename = getImageFileName(date)
	fullPath = imageFolder + str(previousDate.year) + '/' + padzeros(previousDate.month) + '/' + filename
	localpath = "background-image.png"
	print('downloading image ', fullPath)
	with closing(urllib.request.urlopen(fullPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath
	
def padzeros(n):
	"""
	Left pad a number with zeros. 
    """
	return str(n) if n >= 10 else '0'+str(n)
	
def fastSum(start,end):
	refdate = datetime(2017,6,1)
	daysa = (start-refdate).days
	daysb = (end-refdate).days
	days = daysb - daysa + 1
	startsafe = start - timedelta(days = 2)
	endsafe = end
	print(str(daysa) + ", " + str(daysb))
	factor = 1
	dx = None
	dy = None
	while(daysa <= daysb):
		factor *= 2
		if(daysa%factor >= factor/2):
			# add start
			dxx,dyy = loadSimpleFiles(start, start + timedelta(days = factor/2 -1))
			if dx is None:
				dx = dxx
				dy = dyy
			else:
				dx = addbis(dxx,dx)
				dy = addbis(dyy,dy)
			start = start + timedelta(days = factor/2)
			daysa += factor/2
		if(daysb%factor < factor/2 and daysa <= daysb):
			# add end
			dxx,dyy = loadSimpleFiles(end - timedelta(days = factor/2 -1), end)
			if dx is None:
				dx = dxx
				dy = dyy
			else:
				dx = addbis(dxx,dx)
				dy = addbis(dyy,dy)		
			end = end - timedelta(days = factor/2)
			daysb -= factor/2
	s = dateString(startsafe)
	e = dateString(endsafe)
	np.savetxt("dx_average_" + s + "-" + e +".csv", dx, delimiter=",")
	np.savetxt("dy_average_" + s + "-" + e +".csv", dy, delimiter=",")
	
	oldimage = endsafe < datetime(2021,3,12)
	imagescale = 1.25 if oldimage else 1
	im = plotMatrixbis(dx/days, dy/days, imagescale, 'osisaf-test.png')
	#downloadImages(fromDate,date)
	title = padzeros(startsafe.day) + ' ' + monthNames[startsafe.month-1] + " " + str(startsafe.year) + ' to ' + padzeros(endsafe.day) + ' ' + monthNames[endsafe.month-1] + " " + str(endsafe.year)
	filename = 'osisaf-average-' + title + '.png'	
	crop(im, title, imagescale, filename)
	return filename
	
def getImageFileName(date):
	previousDate = date - timedelta(days = 2)
	return 'ice_drift_nh_polstere-625_multi-oi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '1200_combo.png'
	
def prepareImage(filename):
	im = Image.open(filename)
	im = im.convert("RGBA")
	width, height = im.size
	print('image size', width, height)
	pixelmatrix = im.load()
	for row in range(height):
		for col in range(width):
			pixel = pixelmatrix[col, row]
			if not iswater(pixel) and not ismidgrey(pixel):					
				pixelmatrix[col, row] = (255,255,255)
			if ismidgrey(pixel) and not ismidgrey(pixelmatrix[max(col-1,0), row]) and not ismidgrey(pixelmatrix[min(col+1,width-1), row]) and not ismidgrey(pixelmatrix[col, max(row-1,0)]) and not ismidgrey(pixelmatrix[col, min(row+1,height-1)]):
				pixelmatrix[col, row] = (255,255,255)
	im.save('osisaf-average.png')
	
def iswater(pixel):
	return pixel[0] == 4 and pixel[1] == 97 and pixel[2] == 152

def ismidgrey(pixel):
	return pixel[0] == 128 and pixel[1] == 128 and pixel[2] == 128
	
def loadSimpleFiles(start, end):
	filenamex,filenamey = getSimpleFilename(start,end)
	print(filenamex)
	dx = np.loadtxt(filenamex, delimiter=",", dtype=str)
	dx  = np.array([i for i in dx]).astype(float)	
	dy = np.loadtxt(filenamey, delimiter=",", dtype=str)
	dy  = np.array([i for i in dy]).astype(float)
	return dx,dy	

def getSimpleFilename(start,end):
	if(start==end):
		return "data/fast/dx_" + dateString(start) + ".csv", "data/fast/dy_" + dateString(start) + ".csv"
	return "data/fast/dx_" + dateString(start) + "-" + dateString(end) + ".csv", "data/fast/dy_" + dateString(start) + "-" + dateString(end) + ".csv"

def dateString(date, separator = ""):
	return str(date.year) + separator + padzeros(date.month) + separator + padzeros(date.day)

def addbis(dxx,dx):
	for row in range(177):
		for col in range(119):
			if isValid(dx[row,col]):
				if not isValid(dxx[row,col]):
					dxx[row,col] = 0
				dxx[row,col] += dx[row,col]
	
	return dxx

def getSumBis(fromDate, step):
	dx,dy = loadSimpleFiles(fromDate, fromDate + timedelta(days = step/2-1))
	dxx,dyy = loadSimpleFiles(fromDate + timedelta(days = step/2), fromDate + timedelta(days = step-1))
	print(fromDate)
	dxx = addbis(dxx,dx)
	dyy = addbis(dyy,dy)
	return (dxx,dyy)
	
def plotMatrixbis(dx,dy,scale,saveFileName):
	filename = 'osisaf-average.png'
	im = Image.open(filename)
	im = im.convert("RGBA")	
	width, height = im.size
	pixelmatrix = im.load()
	na = np.array(im)
	print(width, height)
	for row in range(177):
		for col in range(119):
			dxx = dx[row,col]
			dyy = dy[row,col]
			if not isValid(dxx) or not isValid(dyy):
				continue
			ii = int(-209*scale+11.5*scale*col)
			jj = int(-602*scale+11.5*scale*row)
			if ii >=0 and jj >= 0 and ii < width and jj < height:
				pixelmatrix[ii, jj] = (0,0,0)
				ptA = (ii,jj)
				ptB = (ii+int(dxx*scale), jj-int(dyy*scale))
				#na = cv2.arrowedLine(na, ptA, ptB, (0,0,0), 1)
				arrowedLine(im,ptA,ptB,scale)
	#im = Image.fromarray(na)
	im.save(saveFileName)
	return im

def isValid(cell):
	if (type(cell) is np.float32):
		return True	
	if((not isnan(cell)) and (type(cell) is np.float64)):
		return True
	return False

def arrowedLine(im, ptA, ptB, scale, width=1, color=(0,0,0)):
	"""Draw line from ptA to ptB with arrowhead at ptB"""
	# Get drawing context
	draw = ImageDraw.Draw(im)
	# Draw the line without arrows
	draw.line((ptA,ptB), width=width, fill=color)
	arrowheadLength = 0.7

	# Now work out the arrowhead
	# = it will be a triangle with one vertex at ptB
	# - it will start at 70% of the length of the line
	# - it will extend 8 pixels either side of the line
	x0, y0 = ptA
	x1, y1 = ptB
	# Now we can work out the x,y coordinates of the bottom of the arrowhead triangle
	xb = arrowheadLength*(x1-x0)+x0
	yb = arrowheadLength*(y1-y0)+y0
	arrowheadWidth = 0.414*(1-arrowheadLength)*sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) # 0.414 
	#if arrowheadWidth >= 2:
	#	print(ptA,ptB,arrowheadWidth)

	# Work out the other two vertices of the triangle
	# Check if line is vertical
	if x0==x1:
		vtx0 = (xb-arrowheadWidth, yb)
		vtx1 = (xb+arrowheadWidth, yb)
	# Check if line is horizontal
	elif y0==y1:
		vtx0 = (xb, yb+arrowheadWidth)
		vtx1 = (xb, yb-arrowheadWidth)
	else:
		alpha = atan2(y1-y0,x1-x0)-90*pi/180
		a = arrowheadWidth*cos(alpha)
		b = arrowheadWidth*sin(alpha)
		vtx0 = (xb+a,yb+b) #(ceil(xb+a) if a > 0 else floor(xb+a), ceil(yb+b) if b > 0 else floor(yb+b))
		vtx1 = (xb-a,yb-b) #(ceil(xb-a) if a < 0 else floor(xb-a), ceil(yb-b) if b < 0 else floor(yb-b))

	#draw.point((xb,yb), fill=(255,0,0))    # DEBUG: draw point of base in red - comment out draw.polygon() below if using this line
	#im.save('DEBUG-base.png')              # DEBUG: save

	# Now draw the arrowhead triangle
	draw.polygon([vtx0, vtx1, ptB], fill=color)
	return im

def crop(im, title, scale, saveFileName):	
	if scale == 1:
		im1 = im.crop((68,30,748,710)) # (318,242,734,658) #(280,220,744,684) #(70,30,744,684)
	elif scale == 1.25:
		im1 = im.crop((76,24,926,874))
	width, height = im1.size	
	pixelmatrix = im1.load()
	print(saveFileName)

	for row in range(int(40*scale)):
		for col in range(width):		
			pixelmatrix[col, row] = (255,255,255)
	printimtext = ImageDraw.Draw(im1)
	fontsize=int(30*scale)
	#font = ImageFont.load_default()
    #font = ImageFont.truetype("/usr/share/fonts/truetype/freefont/arialbd.ttf", fontsize)
	font = ImageFont.truetype("arialbd.ttf", fontsize)
	printimtext.text((5,1), title, (0, 0, 0), font=font)
	
	subtitle1 = "Graphic created by Steven D using"
	subtitle2 = "OSISAF numeric data & background image"
	subtitlefont = ImageFont.truetype("arial.ttf", int(10*scale))
	printimtext.text((int(480*scale),int(5*scale)), subtitle1, (0, 0, 0), font=subtitlefont)
	printimtext.text((int(480*scale),int(20*scale)), subtitle2, (0, 0, 0), font=subtitlefont)

	im1.save(saveFileName)
