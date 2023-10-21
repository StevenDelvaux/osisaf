import numpy as np
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
from dateutil.rrule import rrule, DAILY
import glob
import csv
import sys
import os
import requests
import urllib3
import shutil
import urllib.request
from contextlib import closing
from math import sqrt, sin, cos, pi, floor, isnan, nan, atan, atan2, ceil
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from PIL import Image, ImageGrab, ImageDraw, ImageFont
import cv2
import matplotlib.ticker as ticker
import dropbox
from time import sleep

auto = False
fast = True
updateImage = True
framExport = True
highRes = True
plotAverage = True
plotdays=[40]#[10,30,-1]
enddelta = 0
plotStartYear = 2023
plotEndYear = 2023
framDefinition = 'ateam' #'ateam' '80n' 'flat' 'fjl'
monthLengths = [31,28,31,30,31,30,31,31,30,31,30,31]
monthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
ftpFolder = 'ftp://osisaf.met.no/archive/ice/drift_lr/merged/'
imageFolder = 'https://osisaf.met.no/quicklooks/prod/ice/'

def loadFile(filename):
	#print('loading file', filename)
	f = Dataset(filename, 'r', format="NETCDF4")		
	dx = f.variables['dX'][:]
	dy = f.variables['dY'][:]
	f.close()
	return dx,dy

def loadFileExtended(filename):
	#print('loading file', filename)
	f = Dataset(filename, 'r', format="NETCDF4")		
	dx = f.variables['dX'][:]
	dy = f.variables['dY'][:]
	lat = f.variables['lat'][:]
	lon = f.variables['lon'][:]	
	f.close()
	return dx,dy,lat,lon
	
def fram(filename):
	print('processing',filename)
	dx,dy = loadFile(filename)
	print(dx.shape, dy.shape)
	export = 0
	startIndex = 0
	endIndex = 7 if framDefinition == '80n' else 12 if framDefinition == 'flat' else 13 if framDefinition == 'fjl' else 8
	for k in range(startIndex,endIndex): #range(-1,11):#range(-1,12)#range(-2,5):
		row = (103-1+k) if framDefinition == '80n' else (103-1) if framDefinition == 'flat' else (104-1-k) if framDefinition == 'fjl' else round(102-1+k/3) #103-1#103-1-k#105-1+k
		col = (75-1-k) if framDefinition == '80n' else (65-1+k) if framDefinition == 'flat' else (64-1+k) if framDefinition == 'fjl' else (75-1-k) #65+k#65-1+k#73-1-k
		dxx = dx[0,row,col]
		dyy = dy[0,row,col]
		if not type(dxx) is np.float32 or not type(dyy) is np.float32:
			continue
		export += -dyy if framDefinition == 'flat' else (dxx/3 - dyy) if framDefinition == 'ateam' else (dxx - dyy)
	return export*62.5/2.0

def getFram():
	aframall = []

	for year in range(2023,2024):
		afram = []
		fileList = glob.glob('./data/lr/' + str(year) + '/**/ice_drift_nh_polstere-625_multi-oi*.nc')
		for filename in fileList :
			export = fram(filename)
			#print('export',export)
			afram.append(export)
		aframall.append(afram)
	np.savetxt("fram-osisaf-temp.csv", aframall, delimiter=",", fmt="%.3f")
	#export = fram('./data/lr/2023/01/ice_drift_nh_polstere-625_multi-oi_202301051200-202301071200.nc')
	#print('export',export)
	#afram.append(export)

def getDriftTrack():
	adrift = []
	xx = 0#8*62.5#0
	yy = 0#-8*62.5#0
	counter = 0
	adrift.append([xx,yy,90,0,0,0,93,61,90,0])		
	for year in range(2022,2024):

		fileList = glob.glob('./data/lr/' + str(year) + '/**/ice_drift_nh_polstere-625_multi-oi*.nc')
		for filename in fileList :
			if filename >= './data/lr/2022\\08\\ice_drift_nh_polstere-625_multi-oi_202208271200-202208291200.nc':
				counter += 1
				#if counter > 10:
				#	continue
				print('processing',filename)
				dx,dy,lat,lon = loadFileExtended(filename)
				row = 93 - 1 - int(round(yy/62.5, 0))
				col = 61 - 1 + int(round(xx/62.5, 0))
				dxx = dx[0,row,col]
				dyy = dy[0,row,col]
				if type(dxx) is np.float32 and type(dyy) is np.float32:
					xx += dxx/2
					yy += dyy/2
				excelLatitude = lat[row,col]
				excelLongitude = lon[row,col]
				latitude = 90 - 0.009148928 * sqrt(xx*xx + yy*yy)
				longitude = 45 + atan(yy/xx)/pi*180 if xx > 0 else 225 + atan(yy/xx)/pi*180 if xx < 0 else 45 if yy >= 0 else 315

				print(row+1,col+1, xx, yy, dxx, dyy, row, col, latitude, longitude)			
				adrift.append([xx,yy,latitude,longitude,dxx,dyy,row+1,col+1,excelLatitude,excelLongitude])

		#adriftall.append(adrift)
	np.savetxt("osisaf-drift-temp.csv", adrift, delimiter=",", fmt="%.3f")
	return adrift

def getHighResLandMask():
	landmask = np.genfromtxt('mask125.csv', delimiter=',')	
	print('nsidc land mask',landmask.shape)
	n = landmask.shape[0]  # n is assumed to be an odd number
	offset = 300
	offsetbis = 250
	landmask = 1-landmask[-1:0:-1,:]/2.0
	print('nsidc land mask bis',landmask.shape)
	#exit()
	landmask = landmask[offset:-offsetbis, offset:-offsetbis]
	global landmaskcenter
	landmaskcenter = (n-1)/2 - offset
	print('nsidc land mask cropped',n,landmask.shape,landmaskcenter)
	return landmask

def getNsidcLandMask():
	landmask = np.genfromtxt('landmask_nsidc.csv', delimiter=',')	
	print('nsidc land mask',landmask.shape)
	n = landmask.shape[0]  # n is assumed to be an odd number
	offset = 110
	landmask = landmask[offset:-60, offset:-60]
	global landmaskcenter
	landmaskcenter = (n-1)/2 - offset
	print('nsidc land mask cropped',n,landmask.shape,landmaskcenter)
	return landmask

def insertDataInNsidcMask(landmask, data, dummyvalue = 1, csv = False):		
	print(len(data))
	counter = 0 
	for k in range(len(data)):
		counter += 1
		if counter > 159:
			continue
		coord = data[k]
		print('coord', coord)
		if csv:
			line = coord.split(",")	
			coord = np.array(line[0:2]).astype(float)
			latitude = coord[0]
			longitude = coord[1]
			print('coord bis', coord)
		else:
			xx = coord[0]
			yy = coord[1]
			latitude = 90 - 0.009148928 * sqrt(xx*xx + yy*yy)
			longitude = 45 + atan(yy/xx)/pi*180 if xx > 0 else 225 + atan(yy/xx)/pi*180 if xx < 0 else 45 if yy >= 0 else 315
		rad = (360 if not highRes else 720)*sqrt(2)*sin(pi*(90-latitude)/360) #if not highRes else 
		y = int(round(landmaskcenter+rad*sin(pi*longitude/180.0)))
		x = int(round(landmaskcenter+rad*cos(pi*longitude/180.0)))
		print(k,x,y,latitude,longitude)		
		
		if(x < 0 or y < 0 or x >= landmaskSize or y >= landmaskSize or landmask[x,y] == 0): #or mask[i,j] == -1
			continue				
		#if not type(v) is np.float64:
			#if isnan(refmask[i,j]):
			#	continue
			#v = 0
			#continue
		#if(landmask[x,y] == dummyvalue):
		#	landmask[x,y] = 0
		landmask[x,y] = 5 if csv else 10

	#for x in range(0,landmask.shape[0]):
	#	for y in range(0,landmask.shape[1]):
	#		if(landmask[x,y] == 0 or landmask[x,y] == dummyvalue):
	#			continue
	#		n = floor(landmask[x,y]/counter)
	#		landmask[x,y] = (landmask[x,y] - counter*n)/n
			
	return landmask
	
def plotThickness(landmask,plotTitle,filename,dropboxFilename):
	cdict = {'red':   ((0.0,  0.5, 0.5),
					   (0.001, 0.5, 0.0),
					   (0.5, 0.0, 0.0),
		           	   (1.0,  1, 1)),

         'green':      ((0.0,  0.5, 0.5),
         	           (0.001, 0.5, 0.0),
         	           (0.5, 1.0, 0.0),
         	           (1.0,  0.1, 1)),

         'blue':       ((0.0,  0.5, 0.5),
         	           (0.001, 0.4, 0.4),
					   (0.5, 0.2, 0.2),
         	           (1.0,  0, 1))}
	kleur = LinearSegmentedColormap('BlueRed1', cdict)
	plt.register_cmap(cmap=kleur)

	mask = landmask[100:-60,15:-115]
	print(mask.shape)
	n = landmask.shape[0]
	#plot the relevant map:
	fig2 = plt.imshow(mask, extent=(0,n,0,n), vmin= 0, vmax=11,
           interpolation='nearest', cmap=kleur)
	
	plt.title(plotTitle)
	plt.xticks([])
	plt.yticks([])
	#plt.show()
	plt.savefig(filename)
	#if putPiomasOnDropbox:
	#	plt.savefig(dropboxFilename + '.png')
	#	uploadToDropbox(dropboxFilename + '.png')

def padzeros(n):
	"""
	Left pad a number with zeros. 
    """
	return str(n) if n >= 10 else '0'+str(n)

def getFileName(date):
	previousDate = date - timedelta(days = 2)
	return 'ice_drift_nh_polstere-625_multi-oi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '1200.nc'

def getImageFileName(date):
	previousDate = date - timedelta(days = 2)
	return 'ice_drift_nh_polstere-625_multi-oi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '1200_combo.png'

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

def dateString(date):
	return str(date.year) + padzeros(date.month) + padzeros(date.day)

def getLatestFileInFolder(year):
	listOfFiles = glob.glob('./data/lr/' + str(year) + '/**/ice_drift_nh_polstere-625_multi-oi*.nc')
	return max(listOfFiles, key=os.path.getctime)

def getLatestDateInFolder(year):
	latestFile = getLatestFileInFolder(year)
	print('latest file', latestFile)
	latestDate = datetime(int(latestFile[53:57]), int(latestFile[57:59]), int(latestFile[59:61]))	
	print('here latest', latestDate.day, latestDate.month, latestDate.year)
	return latestDate + timedelta(days = 2)
	
def floatToString(n):
	return "{:.3f}".format(n)
	
def downloadImage(date):
	previousDate = date - timedelta(days = 2)
	filename = getImageFileName(date)
	fullPath = imageFolder + str(previousDate.year) + '/' + padzeros(previousDate.month) + '/' + filename
	localpath = './images/' + filename
	print('downloading image ', fullPath)
	with closing(urllib.request.urlopen(fullPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath	
	
def download(date):
	"""
	Download osisaf ftp file. 
    """
	filename = getFileName(date)
	fullFtpPath = ftpFolder + str(date.year) + '/' + padzeros(date.month) + '/' + filename
	localfolder = './data/lr/' + str(date.year) + '/' + padzeros(date.month)
	localpath = localfolder + '/' + filename
	print('downloading file ', fullFtpPath)
	if not os.path.isdir(localfolder):
		os.makedirs(localfolder, exist_ok=True)
	#r = requests.head(fullFtpPath)
	#fileExists = r.status_code == requests.codes.ok
	#if not fileExists:
	#	return
	#print(fileExists, r.status_code, requests.codes.ok)
	with closing(urllib.request.urlopen(fullFtpPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath	
	
def downloadNewFiles():
	framData = []
	latestDate = getLatestDateInFolder(2023)
	yesterday = datetime.today() - timedelta(days = 1)
	#latestFtpDate = getLatestDateInFolder(ftpFolder)
	date = latestDate + timedelta(days = 1)
	#outFile = open('osisaf-fram-test.csv', 'a', newline='')
	#csvFile = csv.writer(outFile)
	while date < yesterday:
		print('downloading', date, yesterday)
		filename = ''
		try:
			filename = download(date)
		except:
			print('File not found: ', date)
			break
		date = date + timedelta(days = 1)
		framData.append(floatToString(fram(filename)))
		#csvFile.writerow(fram(filename))
	#outFile.close()
	print(framData)
	return framData

def appendToCsvFile(framData):
	if len(framData) == 0:
		return
	with open("osisaf-fram-daily.csv", "a") as myfile:
		myfile.write( ',' + ','.join(framData))

def prepareImage(filename):
	#filename = 'ice_drift_nh_polstere-625_multi-oi_202304111200-202304131200_combo.png' # 'ice_drift_nh_polstere-625_multi-oi_202303071200-202303091200_combo.png'
	im = Image.open(filename)
	im = im.convert("RGBA")
	width, height = im.size
	print('image size', width, height)
	pixelmatrix = im.load()
	for row in range(height):
		for col in range(width):
			pixel = pixelmatrix[col, row]
			#if pixel[0] == 128 and pixel[1] == 0 and pixel[2] == 128:
			#if pixel[0] == 128 and pixel[1] == 0 and pixel[2] == 128:				
			#pixelmatrix[col, row] = (255,255,255)
			#if pixel[0] != 128 and pixel[0] != 255 and pixel[0] == pixel[1] and pixel[0] == pixel[2]:
			#if (pixel[0] != 128 or pixel[1] != 128) and pixel[0] == pixel[2]:
			if not iswater(pixel) and not ismidgrey(pixel):					
				pixelmatrix[col, row] = (255,255,255)
			if ismidgrey(pixel) and not ismidgrey(pixelmatrix[max(col-1,0), row]) and not ismidgrey(pixelmatrix[min(col+1,width-1), row]) and not ismidgrey(pixelmatrix[col, max(row-1,0)]) and not ismidgrey(pixelmatrix[col, min(row+1,height-1)]):
				pixelmatrix[col, row] = (255,255,255)
	im.save('osisaf-average.png')	

def iswater(pixel):
	return pixel[0] == 4 and pixel[1] == 97 and pixel[2] == 152

def ismidgrey(pixel):
	return pixel[0] == 128 and pixel[1] == 128 and pixel[2] == 128
	
def downloadImages(fromDate, toDate):
	date = fromDate
	while date <= toDate:
		downloadImage(date)
		date = date + timedelta(days = 1)

def cropImages(fromDate, toDate):
	date = fromDate
	while date <= toDate:
		cropImage(date)
		date = date + timedelta(days = 1)
	
def plotImage(date):
	filename = getFileName(date)
	localpath = './data/lr/' + str(date.year) + '/' + padzeros(date.month) + '/' + filename
	print('gonna plot image',filename)
	dx,dy = loadFile(localpath)
	print(dx.shape, dy.shape)
	plotMatrix(dx,dy,'osisaf-test.png')

def plotMatrix(dx,dy,saveFileName):
	filename = 'osisaf-average.png'
	im = Image.open(filename)
	im = im.convert("RGBA")	
	width, height = im.size
	pixelmatrix = im.load()
	na = np.array(im)
	print(width, height)
	
	for row in range(177):
		for col in range(119):
			dxx = dx[0,row,col]
			dyy = dy[0,row,col]
			if not type(dxx) is np.float32 or not type(dyy) is np.float32:
				continue
			ii = int(-209+11.5*col)
			jj = int(-602+11.5*row)
			if ii >=0 and jj >= 0 and ii < width and jj < height:
				pixelmatrix[ii, jj] = (0,0,0)
				scale = 1.0
				ptA = (ii,jj)
				ptB = (ii+int(dxx/scale), jj-int(dyy/scale))
				#na = cv2.arrowedLine(na, ptA, ptB, (0,0,0), 1)
				arrowedLine(im,ptA,ptB)
	#im = Image.fromarray(na)
	im.save(saveFileName)
	return im

def plotMatrixbis(dx,dy,saveFileName):
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
			ii = int(-209+11.5*col)
			jj = int(-602+11.5*row)
			if ii >=0 and jj >= 0 and ii < width and jj < height:
				pixelmatrix[ii, jj] = (0,0,0)
				scale = 1.0
				ptA = (ii,jj)
				ptB = (ii+int(dxx/scale), jj-int(dyy/scale))
				#na = cv2.arrowedLine(na, ptA, ptB, (0,0,0), 1)
				arrowedLine(im,ptA,ptB)
	#im = Image.fromarray(na)
	im.save(saveFileName)
	return im
	
def arrowedLine(im, ptA, ptB, width=1, color=(0,0,0)):
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
	
def cropImage(date):
	previousDate = date - timedelta(days = 2)
	filename = getImageFileName(date)
	fullPath = './images/' + filename
	im = Image.open(fullPath)
	im = im.convert("RGBA")
	title = padzeros(previousDate.day) + ' ' + monthNames[previousDate.month-1] + " " + str(previousDate.year) + ' to ' + padzeros(date.day) + ' ' + monthNames[date.month-1] + " " + str(date.year)	
	saveFileName = './imagecrop/' + filename
	crop(im, title, saveFileName)

def crop(im, title, saveFileName):	
	im1 = im.crop((68,30,748,710)) # (318,242,734,658) #(280,220,744,684) #(70,30,744,684)
	width, height = im1.size	
	pixelmatrix = im1.load()
	print(saveFileName)

	for row in range(40):
		for col in range(width):		
			pixelmatrix[col, row] = (255,255,255)
	printimtext = ImageDraw.Draw(im1)
	fontsize=30 
	font = ImageFont.truetype("arialbd.ttf", fontsize)
	printimtext.text((5,1), title, (0, 0, 0), font=font)

	im1.save(saveFileName)

def isValid(cell):
	if (type(cell) is np.float32):
		return True	
	if((not isnan(cell)) and (type(cell) is np.float64)):
		return True
	return False
	
def add(dxx,dx):
	for row in range(177):
		for col in range(119):
			if isValid(dx[0,row,col]):
				if not isValid(dxx[0,row,col]):
					dxx[0,row,col] = 0
				dxx[0,row,col] += dx[0,row,col]
	
	return dxx

def addbis(dxx,dx):
	for row in range(177):
		for col in range(119):
			if isValid(dx[row,col]):
				if not isValid(dxx[row,col]):
					dxx[row,col] = 0
				dxx[row,col] += dx[row,col]
	
	return dxx

def getSum(fromDate, toDate):
	date = fromDate
	filename = getFileName(date)
	localpath = './data/lr/' + str(date.year) + '/' + padzeros(date.month) + '/' + filename
	dxx,dyy = loadFile(localpath)
	date = date + timedelta(days = 1)
	
	while date <= toDate:
		print(date)
		filename = getFileName(date)
		localpath = './data/lr/' + str(date.year) + '/' + padzeros(date.month) + '/' + filename
		dx,dy = loadFile(localpath)
		dxx = add(dxx,dx)
		dyy = add(dyy,dy)
		date = date + timedelta(days = 1)
	return (dxx,dyy)
	
def getSumBis(fromDate, step):
	dx,dy = loadSimpleFiles(fromDate, fromDate + timedelta(days = step/2-1))
	dxx,dyy = loadSimpleFiles(fromDate + timedelta(days = step/2), fromDate + timedelta(days = step-1))
	print(fromDate)
	dxx = addbis(dxx,dx)
	dyy = addbis(dyy,dy)
	return (dxx,dyy)
	
def plotPlot(delta, startyear, endyear, enddelta):
	print('inside plot', delta, enddelta, startyear, endyear)
	#prepareImage()
	#exit()
	years = endyear - startyear + 1
	dxx=None
	dyy=None
	enddate = datetime.today() - timedelta(days = 1)
	for year in range(startyear,endyear+1):
		date = datetime(year, 3, 20)	
		fromDate = datetime(year, 3, 11)
		if delta > 0:
			date = datetime(year, enddate.month, enddate.day)
			fromDate = date - timedelta(days = delta-2)
			date = date - timedelta(days = enddelta)
		elif delta == -1:
			date = datetime(year, enddate.month, enddate.day)
			fromDate = datetime(year-1, 10, 3)		
		days = (date - fromDate).days + 1
		dxxx,dyyy = getSum(fromDate, date)
		if dxx is None:
			dxx = dxxx
			dyy = dyyy
		else:
			dxx = add(dxx,dxxx)
			dyy = add(dyy,dyyy)
			
	im = plotMatrix(dxx/days/years, dyy/days/years, 'osisaf-test.png')
	#downloadImages(fromDate,date)
	fromDateBefore = fromDate - timedelta(days = 2)
	dateBefore = date #- timedelta(days = 1)
	title = str(fromDateBefore.day) + ' ' + monthNames[fromDateBefore.month-1] + " " + str(startyear-1 if delta == -1 else startyear) + ' to ' + str(dateBefore.day) + ' ' + monthNames[dateBefore.month-1] + " " + str(dateBefore.year)	
	aux = '' if auto else ('-' + str(endyear)) if startyear == endyear else ('-' + str(startyear) + '-' + str(endyear))
	crop(im, title, 'osisaf-average-' + ('last-' if enddelta == 0 else '') + str(delta) + ('-days' if enddelta == 0 else '-to-' + str(enddelta) + '-days') + aux + '.png' if delta > 0 else 'osisaf-average' + aux + '-since-october.png' if delta == -1 else 'osisaf-avg.png')
	#cropImages(fromDate,date)	
	#plotImage(date)

def plotCumSum(ax, lines, dates, idx, label, color, memoryData, linewidth=1):
	line = lines[idx].replace(',,', ',0,').replace(',,', ',0,')
	line = line.split(",")
	if len(memoryData) > 0:	
		merged = np.append(memoryData,line[1:274])
	
		row =  np.array([i.lstrip() for i in np.array(merged)])
		cum = np.cumsum(row.astype(float))/1000000.0
		ax.plot(dates, cum, label=label, color=color, linewidth=linewidth);	
	return line[274:366]
	
def plotFramGraph(inputFileName, outputFileName, title, ymax):
	fig, ax = plt.subplots(figsize=(8, 5))
	dates = np.arange(1,366)	
	with open(inputFileName, 'r') as f:
		lines = f.readlines()
	memoryData = []
	memoryData = plotCumSum(ax, lines, dates, -9, '2014/15', (0.0,0.69,0.94), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -8, '2015/16', (0,0.69,0.31), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -7, '2016/17', (0.57,0.82,0.31), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -6, '2017/18', (1.0,0.75,0), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -5, '2018/19', (0.9,0.4,0.05), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -4, '2019/20', (1.0,0.5,0.5), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -3, '2020/21', (0.58,0.54,0.33), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -2, '2021/22', (0.4,0,0.2), memoryData)
	#plotCumSum(ax, lines, dates, -1, '2022/23', (1.0,0,0), linewidth=2)
	
	#line = lines[-1].split(",")
	#line21 = np.array([i.lstrip() for i in np.array(line[1:366])])
	#line21[364] = line21[364].replace("\r2022", "" )	
	#cum21 = np.cumsum(line21.astype(float))
	#ax.plot(dates, cum21, label="2021", color=(0.4,0,0.2));
	
	line = lines[-1].split(",")
	line22 = np.append(memoryData, (np.array([i.lstrip() for i in np.array(line[1:])])))
	numberOfDays = len(line22)
	print(numberOfDays)
	cum22 = np.cumsum(line22.astype(float))/1000000.0
	padded22 = np.pad(cum22, (0, 365 - numberOfDays), 'constant', constant_values=(np.nan,))
	ax.plot(dates, padded22, label="2022/23", color=(1.0,0,0), linewidth=3);
	
	#ax.set_xlabel("day")
	ax.set_ylabel("10$^6$ km$^2$")
	ax.set_title(title)
	plt.text(135,1.1,'calculated from OSI-405')
	ax.legend(loc=2, prop={'size': 8})
	ax.axis([0, 365, 0, ymax])
	ax.grid(True);
	
	months = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep']
	ax.set_xticks([0,31,61,92,123,151,182,212,243,273,304,335,365], ['', '', '', '', '', '', '', '', '', '', '', '', ''])
	ax.xaxis.set_minor_locator(ticker.FixedLocator([15.5,46,76.5,107.5,137,166.5,197,227.5,258,288.5,319.5,350]))
	ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months))
	ax.tick_params(which='minor', length=0)
		
	fig.savefig(outputFileName)

def uploadToDropbox(filename):	
	dropbox_access_token= "Xdyuu7yOogcAAAAAAAAAAfDLocTCuWt_byjRA3lmVgpjLZScP7S9jHwrvGRCgl-6"
	client = dropbox.Dropbox(dropbox_access_token)	
	print("[SUCCESS] dropbox account linked")
	
	dropbox_path= "/" + filename
	computer_path= "./" + filename
	client.files_upload(open(computer_path, "rb").read(), dropbox_path, mode=dropbox.files.WriteMode.overwrite)
	print("[UPLOADED] {}".format(computer_path))

def addNans(dx):
	for i in range(dx.shape[0]):
		for j in range(dx.shape[1]):
			if not isValid(dx[i,j]):
				dx[i,j] = nan
				
def generatePrecalculated(startdate, enddate, step): # step is a power of 2
	counter = 0
	while startdate + timedelta(days = step-1) <= enddate:	
		counter += 1
		#if counter > 2:
		#	exit()
		if step==1:
			print(startdate)
			dxxx,dyyy = getSum(startdate, startdate)
			dxxx = dxxx[0,:,:]
			dyyy = dyyy[0,:,:]
		elif step==2:
			dxxx,dyyy = getSum(startdate, startdate + timedelta(days = 1))
			dxxx = dxxx[0,:,:]
			dyyy = dyyy[0,:,:]
		else:
			dxxx,dyyy = getSumBis(startdate, step)	

		#print(dxxx.shape)
		filex,filey = getSimpleFilename(startdate, startdate + timedelta(days = step-1))
		#print(dxxx[0,0])
		#print(type(dxxx[0,0]))
		#print(dxxx[89,26])
		#print(type(dxxx[89,26]))
		addNans(dxxx)
		addNans(dyyy)
		np.savetxt(filex, dxxx, delimiter=",")
		np.savetxt(filey, dyyy, delimiter=",")
		startdate = startdate + timedelta(days = step)

def fastSum(start,end):
	refdate = datetime(2017,6,1)
	daysa = (start-refdate).days
	daysb = (end-refdate).days
	days = daysb - daysa + 1
	startsafe = start
	endsafe = end
	print(str(daysa) + ", " + str(daysb))
	factor = 1
	dx = None
	dy = None
	while(daysa < daysb):
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
		if(daysb%factor < factor/2 and daysa < daysb):
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
	
	im = plotMatrixbis(dx/days, dy/days, 'osisaf-test.png')
	#downloadImages(fromDate,date)
	title = s + "-" + e	
	crop(im, title, 'osisaf-average-' + title + '.png')

if auto:
	plotdays = [10,30,-1]
	framData = downloadNewFiles()
	tomorrow = datetime.today() + timedelta(days = 1)
	if updateImage:
		yesterday = datetime.today() - timedelta(days = 1)
		filename = downloadImage(yesterday)
		prepareImage(filename)
	
	for delta in plotdays:
		plotPlot(delta, tomorrow.year, tomorrow.year, 0)	
	appendToCsvFile(framData)	
elif(fast):	
	step = 32
	startdate = datetime(2017,6,1)
	yesterday = datetime.today() - timedelta(days = 1)
	enddate = datetime(yesterday.year, yesterday.month, yesterday.day)
	
	fastSum(enddate - timedelta(days = 30),enddate)
	exit()
	
	generatePrecalculated(startdate, enddate, step)
	exit()
	
	counter = 0
	while startdate < enddate:
		counter += 1
		if counter > 4:
			exit()
		#dxxx,dyyy = getSum(startdate, startdate + timedelta(days = 1))
		#dxxx = dxxx[0,:,:]
		#dyyy = dyyy[0,:,:]
		#print(dxxx.shape)
		#filex,filey = getSimpleFilename(startdate, startdate + timedelta(days = 1))
		#print(dxxx[0,0])
		#print(type(dxxx[0,0]))
		#print(dxxx[89,26])
		#print(type(dxxx[89,26]))
		#addNans(dxxx)
		#addNans(dyyy)
		#np.savetxt(filex, dxxx, delimiter=",")
		#np.savetxt(filey, dyyy, delimiter=",")
		#startdate = startdate + timedelta(days = 2)
		#sleep(5)
		#print("loading")
		#data = np.loadtxt(filex, delimiter=",", dtype=str)
		#data = np.array([i for i in data]).astype(float)
		#print(type(data))
		#print(data.shape)
		#print(data[0,0])
		#print(isnan(data[0,0]))
		#print(type(data[0,0]))
		#print(data[89,26])
		#print(isnan(data[89,26]))
		#print(type(data[89,26]))
	#exit()
elif plotAverage:
	for delta in plotdays:
		plotPlot(delta, plotStartYear, plotEndYear, enddelta)
	exit()
	
if framExport:
	getFram()
	
	if auto:
		dropboxFileName = "osisaf-fram-daily"
		dropboxFileNameCsv = dropboxFileName + ".csv"
		dropboxFileNamePng = dropboxFileName + ".png"
		#plotFramGraph(dropboxFileNameCsv, dropboxFileNamePng, "OSISAF cumulative sea ice extent export through Fram Strait", 1.15)

		uploadToDropbox(dropboxFileNameCsv)
		uploadToDropbox(dropboxFileNamePng)
	
	exit()
	
driftPath = getDriftTrack()

landmask = getHighResLandMask() if highRes else getNsidcLandMask()
dummyvalue = 1
landmask = landmask*dummyvalue
landmaskSize = landmask.shape[0]
	
with open('BuoyData.csv', 'r') as f:
	lines = f.readlines()
	insertDataInNsidcMask(landmask,lines,csv=True)
landmask = insertDataInNsidcMask(landmask,driftPath)

plotThickness(landmask,'OSISAF simulated drift path 27 Aug 2022 to 27 Jan 2023','osisaf-drift-path.png','')
