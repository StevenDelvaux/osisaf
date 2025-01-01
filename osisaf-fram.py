import numpy as np
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
import glob
import os
import requests
import shutil
import urllib.request
from contextlib import closing
from math import sqrt, sin, cos, pi, floor, isnan, nan, atan, atan2, ceil
import matplotlib.pyplot as plt
from PIL import Image, ImageGrab, ImageDraw, ImageFont
import matplotlib.ticker as ticker
import dropbox
from time import sleep
from decouple import config

auto = True
updateImage = True
framExport = True
plotdays=[10,30]
enddelta = 0
framDefinition = 'ateam' #'ateam' '80n' 'flat' 'fjl'
monthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
ftpFolder = 'ftp://osisaf.met.no/archive/ice/drift_lr/merged/'
imageFolder = 'https://osisaf.met.no/quicklooks/prod/ice/'
putOnDropbox = True

def loadFile(filename):
	#print('loading file', filename)
	f = Dataset(filename, 'r', format="NETCDF4")		
	dx = f.variables['dX'][:]
	dy = f.variables['dY'][:]
	f.close()
	return dx,dy
	
def fram(filename, framDefinition):
	dx,dy = loadFile(filename)
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
	
def padzeros(n):
	"""
	Left pad a number with zeros. 
    """
	return str(n) if n >= 10 else '0'+str(n)

def getFileName(date, north=True):
	previousDate = date - timedelta(days = 2)
	return 'ice_drift_' + ('nh' if north else 'sh') + '_polstere-625_multi-oi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '1200.nc'

def getImageFileName(date, north=True):
	previousDate = date - timedelta(days = 2)
	if north:
		return 'ice_drift_nh_polstere-625_multi-oi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200-' + str(date.year) + padzeros(date.month) + padzeros(date.day) + '1200_combo.png'
	else:
		return 'ice_conc_sh_polstere-100_multi_' + str(previousDate.year) + padzeros(previousDate.month) + padzeros(previousDate.day) + '1200_qlook.png'	

def loadSimpleFiles(start, end, north = True):
	filenamex,filenamey = getSimpleFilename(start, end, north)
	print(filenamex)
	dx = np.loadtxt(filenamex, delimiter=",", dtype=str)
	dx  = np.array([i for i in dx]).astype(float)	
	dy = np.loadtxt(filenamey, delimiter=",", dtype=str)
	dy  = np.array([i for i in dy]).astype(float)
	return dx,dy	

def getSimpleFilename(start, end, north):
	prefix = './data/' + ('fast' if north else 'antarctic') + '/'
	if(start==end):
		suffix = dateString(start) + ".csv"
	else:
		suffix = dateString(start) + "-" + dateString(end) + ".csv"
	return prefix + "dx_" + suffix, prefix + "dy_" + suffix

def dateString(date):
	return str(date.year) + padzeros(date.month) + padzeros(date.day)

def getLatestFileInFolder():
	listOfFiles = glob.glob('./data/antarctic/*.csv')
	return max(listOfFiles)

def getLatestDateInFolder():
	latestFile = getLatestFileInFolder()
	print('latest file',  latestFile)
	latestDate = datetime(int(latestFile[-12:-8]), int(latestFile[-8:-6]), int(latestFile[-6:-4]))	
	print('here latest', latestDate.day, latestDate.month, latestDate.year)
	#return datetime(2024,10,12) # todo temp
	return latestDate
	
def floatToString(n):
	return "{:.3f}".format(n)
	
def downloadImage(date, north=True):
	previousDate = date - timedelta(days = 2)
	filename = getImageFileName(date, north)
	fullPath = imageFolder + str(previousDate.year) + '/' + padzeros(previousDate.month) + '/' + filename
	localpath = './background-image-' + ('arctic' if north else 'antarctic') + '.png'
	print('downloading image ', fullPath)
	with closing(urllib.request.urlopen(fullPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath
	
def download(date, north=True):
	"""
	Download osisaf ftp file. 
    """
	filename = getFileName(date, north)
	fullFtpPath = ftpFolder + str(date.year) + '/' + padzeros(date.month) + '/' + filename
	localfolder = './data/lr' 
	localpath = localfolder + '/' + filename
	print('downloading file ', fullFtpPath)
	if not os.path.isdir(localfolder):
		os.makedirs(localfolder, exist_ok=True)
	with closing(urllib.request.urlopen(fullFtpPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	return localpath	
	
def downloadNewFiles():
	framData = []
	framDataFjl = []
	latestDate = getLatestDateInFolder()
	yesterday = datetime.today() - timedelta(days = 1)
	date = latestDate + timedelta(days = 1)
	while date < yesterday:
		print('downloading', date, yesterday)
		filename = ''
		try:
			filename = download(date)
			filenameAntarctic = download(date, False)
		except:
			print('File not found: ', date)
			break
		date = date + timedelta(days = 1)
		framData.append(floatToString(fram(filename,'ateam')))
		framDataFjl.append(floatToString(fram(filename,'fjl')))
	print(framData)
	return framData, framDataFjl

def appendToCsvFile(framData, filename):
	if len(framData) == 0:
		return
	with open(filename, "a") as myfile:
		myfile.write( ',' + ','.join(framData))

def prepareImage(filename, filenameToSave):
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
	im.save(filenameToSave)	

def iswater(pixel):
	return pixel[0] == 4 and pixel[1] == 97 and pixel[2] == 152

def ismidgrey(pixel):
	return pixel[0] == 128 and pixel[1] == 128 and pixel[2] == 128
	
def downloadImages(fromDate, toDate):
	date = fromDate
	while date <= toDate:
		downloadImage(date)
		date = date + timedelta(days = 1)

def plotMatrixbis(dx, dy, filename, saveFileName):
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
	
	subtitle1 = "Graphic created by Steven D using"
	subtitle2 = "OSISAF numeric data & background image"
	subtitlefont = ImageFont.truetype("arial.ttf", 10)
	printimtext.text((480,5), subtitle1, (0, 0, 0), font=subtitlefont)
	printimtext.text((480,20), subtitle2, (0, 0, 0), font=subtitlefont)

	im1.save(saveFileName)

def isValid(cell):
	if (type(cell) is np.float32):
		return True	
	if((not isnan(cell)) and (type(cell) is np.float64)):
		return True
	return False
	
def add(dxx, dx, north=True):
	for row in range(177 if north else 131):
		for col in range(119 if north else 125):
			if isValid(dx[0,row,col]):
				if not isValid(dxx[0,row,col]):
					dxx[0,row,col] = 0
				dxx[0,row,col] += dx[0,row,col]
	
	return dxx

def addbis(dxx,dx, north = True):
	for row in range(177 if north else 131):
		for col in range(119 if north else 125):
			if isValid(dx[row,col]):
				if not isValid(dxx[row,col]):
					dxx[row,col] = 0
				dxx[row,col] += dx[row,col]
	
	return dxx

def loadFileForDate(date, north=True):
	filename = getFileName(date, north)
	localpath = './data/lr/' + filename
	dxx,dyy = loadFile(localpath)
	return (dxx,dyy)

def getSum(fromDate, toDate, north=True):
	date = fromDate
	print(date)
	filename = getFileName(date, north)
	localpath = './data/lr/' + filename
	dxx,dyy = loadFile(localpath)
	date = date + timedelta(days = 1)
	
	while date <= toDate:
		print(date)
		filename = getFileName(date, north)
		localpath = './data/lr/' + filename
		dx,dy = loadFile(localpath)
		dxx = add(dxx, dx, north)
		dyy = add(dyy, dy, north)
		date = date + timedelta(days = 1)
	return (dxx,dyy)
	
def getSumBis(fromDate, step):
	dx,dy = loadSimpleFiles(fromDate, fromDate + timedelta(days = step/2-1))
	dxx,dyy = loadSimpleFiles(fromDate + timedelta(days = step/2), fromDate + timedelta(days = step-1))
	print(fromDate)
	dxx = addbis(dxx,dx)
	dyy = addbis(dyy,dy)
	return (dxx,dyy)
	
def plotCumSum(ax, lines, dates, idx, label, color, memoryData, linewidth=1):
	line = lines[idx].replace(',,', ',0,').replace(',,', ',0,')
	line = line.split(",")
	if len(memoryData) > 0:	
		merged = np.append(memoryData,line[1:274]) # start plot on 1 October
	
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
	memoryData = plotCumSum(ax, lines, dates, -11, '2014/15', (0.0,0.69,0.94), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -10, '2015/16', (0,0.69,0.31), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -9, '2016/17', (0.57,0.82,0.31), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -8, '2017/18', (1.0,0.75,0), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -7, '2018/19', (0.9,0.4,0.05), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -6, '2019/20', (1.0,0.5,0.5), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -5, '2020/21', (0.58,0.54,0.33), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -4, '2021/22', (0.4,0,0.2), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -3, '2022/23', (0.0,0.69,0.94), memoryData)
	memoryData = plotCumSum(ax, lines, dates, -2, '2023/24', (0,0.44,0.75), memoryData)
	# numberOfDays = len(memoryData)
	# print(numberOfDays)
	# cum23 = np.cumsum(np.array(memoryData, dtype='float32'))/1000000.0
	# padded23 = np.pad(cum23, (0, 365 - numberOfDays), 'constant', constant_values=(np.nan,))	
	# ax.plot(dates, padded23, label="2024/25", color=(1.0,0,0), linewidth=3);
	
	line = lines[-1].split(",")
	line23 = np.append(memoryData, (np.array([i.lstrip() for i in np.array(line[1:])])))
	numberOfDays = len(line23)
	print(numberOfDays)
	cum23 = np.cumsum(line23.astype(float))/1000000.0
	padded23 = np.pad(cum23, (0, 365 - numberOfDays), 'constant', constant_values=(np.nan,))
	ax.plot(dates, padded23, label="2024/25", color=(1.0,0,0), linewidth=3);
	
	
	ax.set_ylabel("10$^6$ km$^2$")
	ax.set_title(title)
	plt.text(95,ymax-0.05,'Graphic created by Steven D. using OSI-405 numeric data', fontsize=8)	
	ax.legend(loc=2, prop={'size': 8})
	ax.axis([0, 365, 0, ymax])
	ax.grid(True);
	
	months = ['Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep']
	ax.set_xticks([0,31,61,92,123,151,182,212,243,273,304,335,365], ['', '', '', '', '', '', '', '', '', '', '', '', ''])
	ax.xaxis.set_minor_locator(ticker.FixedLocator([15.5,46,76.5,107.5,137,166.5,197,227.5,258,288.5,319.5,350]))
	ax.xaxis.set_minor_formatter(ticker.FixedFormatter(months))
	ax.tick_params(which='minor', length=0)
		
	print('gonna save fig: ' + outputFileName)
	fig.savefig(outputFileName)

def addNans(dx):
	for i in range(dx.shape[0]):
		for j in range(dx.shape[1]):
			if not isValid(dx[i,j]):
				dx[i,j] = nan
				
def generatePrecalculated(startdate, enddate, step, north): # step is a power of 2
	print('generatePrecalculated', startdate, enddate, step)
	print(startdate)
	print(enddate)
	print(step)
	counter = 0
	while startdate + timedelta(days = step-1) <= enddate:
		print('generatePrecalculated inner')
		counter += 1
		if step==1:
			print(startdate)
			dxxx,dyyy = loadFileForDate(startdate, north)
			dxxx = dxxx[0,:,:]
			dyyy = dyyy[0,:,:]
		else:
			dxxx,dyyy = getSumBis(startdate, step)	

		filex,filey = getSimpleFilename(startdate, startdate + timedelta(days = step-1), north)
		addNans(dxxx)
		addNans(dyyy)
		print("saving: " + filex)
		print("saving: " + filey)
		np.savetxt(filex, dxxx, delimiter=",")
		np.savetxt(filey, dyyy, delimiter=",")
		startdate = startdate + timedelta(days = step)

def fastSum(start, end, filename):
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
	#np.savetxt("dx_average_" + s + "-" + e +".csv", dx, delimiter=",")
	#np.savetxt("dy_average_" + s + "-" + e +".csv", dy, delimiter=",")
	
	im = plotMatrixbis(dx/days, dy/days, 'osisaf-average.png', 'osisaf-test.png')
	title = padzeros(startsafe.day) + ' ' + monthNames[startsafe.month-1] + " " + str(startsafe.year) + ' to ' + padzeros(endsafe.day) + ' ' + monthNames[endsafe.month-1] + " " + str(endsafe.year)
	crop(im, title, filename)

def updateFramGraphs(framDef):
	framDefinition = framDef

	dropboxFileName = 'osisaf-fram-' + ('fjl-' if framDefinition == 'fjl' else '') + 'daily'
	dropboxFileNameCsv = dropboxFileName + ".csv"
	dropboxFileNamePng = dropboxFileName + ".png"
	plotFramGraph(dropboxFileNameCsv, dropboxFileNamePng, "OSISAF cumulative sea ice extent export through Fram" + ('-FJL' if framDefinition == 'fjl' else ' Strait'), 1.6 if framDefinition == 'fjl' else 1.15)
	if putOnDropbox:
		uploadToDropbox([dropboxFileNameCsv, dropboxFileNamePng])

def plotMatrixAntarctic(dx, dy, filename, saveFileName):
	im = Image.open(filename)
	im = im.convert("RGBA")	
	width, height = im.size
	pixelmatrix = im.load()
	na = np.array(im)
	print(width, height)
	
	for row in range(131):
		for col in range(125):
			dxx = dx[row,col]
			dyy = dy[row,col]
			
			if not isValid(dxx) or not isValid(dyy):
				continue
			ii = int(607+9.35*(col-62))
			jj = int(745+9.35*(row-68))
			if ii >=0 and jj >= 0 and ii < width and jj < height:
				pixelmatrix[ii, jj] = (0,0,0)
				scale = 1.0
				ptA = (ii,jj)
				ptB = (ii+int(dxx/scale), jj-int(dyy/scale))
				arrowedLine(im,ptA,ptB)
	im.save(saveFileName)
	return im
	
def sumAntarctic(start, end, filename):
	fromDateBefore = start - timedelta(days = 2)
	days = (end - start).days + 1
	date = start
	dx = None
	dy = None	
	while date <= end:
		dxx,dyy = loadSimpleFiles(date, date, False)
		if dx is None:
			dx = dxx
			dy = dyy
		else:
			dx = addbis(dxx, dx, False)
			dy = addbis(dyy, dy, False)
		date = date + timedelta(days = 1)
	
	im = plotMatrixAntarctic(dx/days, dy/days, 'osisaf-antarctic-average.png', 'osisaf-antarctic-test.png')
	title = str(fromDateBefore.day) + ' ' + monthNames[fromDateBefore.month-1] + " " + str(fromDateBefore.year) + ' to ' + str(end.day) + ' ' + monthNames[end.month-1] + " " + str(end.year)	
	cropAntarctic(im, title, filename)
	
def cropAntarctic(im, title, saveFileName):	
	im1 = im.crop((157,163,1133,1203)) #21,60,1191,1330#1392,1357)) #((68,30,748,710)) # (318,242,734,658) #(280,220,744,684) #(70,30,744,684)
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
	
	subtitle1 = "Graphic created by Steven D using"
	subtitle2 = "OSISAF numeric data & background image"
	subtitlefont = ImageFont.truetype("arial.ttf", 10)
	printimtext.text((730,5), subtitle1, (0, 0, 0), font=subtitlefont)
	printimtext.text((730,20), subtitle2, (0, 0, 0), font=subtitlefont)

	im1.save(saveFileName)
	
def prepareImageAntarctic(filename, filenameToSave):
	im = Image.open(filename)
	im = im.convert("RGBA")
	width, height = im.size
	print('image size', width, height)
	pixelmatrix = im.load()
	for row in range(height):
		for col in range(width):
			pixel = pixelmatrix[col, row]
			if col > 1012 and row < 175:			
				pixelmatrix[col, row] = (4,98,153)
			if col > 168 and col < 480 and row > 1193:			
				pixelmatrix[col, row] = (4,98,153)
	im.save(filenameToSave)

def uploadToDropbox(filenames):
	dropbox_access_token = config('DROPBOX_ACCESS_TOKEN')
	app_key = config('APP_KEY')
	app_secret = config('APP_SECRET')
	oauth2_refresh_token = config('OAUTH2_REFRESH_TOKEN')
	client = dropbox.Dropbox(oauth2_access_token=dropbox_access_token,app_key=app_key,app_secret=app_secret,oauth2_refresh_token=oauth2_refresh_token)
	print("[SUCCESS] dropbox account linked")
	
	for computer_path in filenames:
		print("[UPLOADING] {}".format(computer_path))
		dropbox_path= "/" + computer_path
		client.files_upload(open(computer_path, "rb").read(), dropbox_path, mode=dropbox.files.WriteMode.overwrite)
		print("[UPLOADED] {}".format(computer_path))
	
def downloadFromDropbox(filenames):
	dropbox_access_token = config('DROPBOX_ACCESS_TOKEN')
	app_key = config('APP_KEY')
	app_secret = config('APP_SECRET')
	oauth2_refresh_token = config('OAUTH2_REFRESH_TOKEN')
	client = dropbox.Dropbox(oauth2_access_token=dropbox_access_token,app_key=app_key,app_secret=app_secret,oauth2_refresh_token=oauth2_refresh_token)
	print("[SUCCESS] dropbox account linked")
	
	for dropbox_path in filenames:
		print("[DWONLOADING] {}".format(dropbox_path))
		computer_path= dropbox_path
		with open(computer_path, 'wb') as f:
			metadata, res = client.files_download(path= "/" + dropbox_path)
			f.write(res.content)
		print("[DOWNLOADED] {}".format(dropbox_path))
	
if auto:
	plotdays = [10,30]
	latestDate = getLatestDateInFolder()
	framData,framDataFjl = downloadNewFiles()
	tomorrow = datetime.today() + timedelta(days = 1)
	yesterday = datetime.today() - timedelta(days = 1)
	if updateImage:	
		filenameImage = downloadImage(yesterday)
		prepareImage(filenameImage, 'osisaf-average.png')
		os.remove(filenameImage)
		
		filenameImageAntarctic = downloadImage(yesterday, False)
		prepareImageAntarctic(filenameImageAntarctic, 'osisaf-antarctic-average.png')
		os.remove(filenameImageAntarctic)
		
	if framExport:
		framFilename = 'osisaf-fram-daily.csv'
		framFjlFilename = 'osisaf-fram-fjl-daily.csv'
		downloadFromDropbox([framFilename, framFjlFilename])
		appendToCsvFile(framData, framFilename)
		appendToCsvFile(framDataFjl, framFjlFilename)
	print('inside precalculate', yesterday, latestDate)
	refdate = datetime(2017,6,1)
	date = latestDate + timedelta(days = 1) #
	hasnew = False
	while date < yesterday: 
		print(date)
		generatePrecalculated(date, date, 1, False) # generate Antarctic files
		
		days = (date - refdate).days + 1
		print(days)
		factor = 1
		hasnew = True
		while days%factor == 0: # generate Arctic files
			print(factor)
			sleep(5)
			generatePrecalculated(date - timedelta(days = factor - 1), date, factor, True)
			factor *= 2			
		date = date + timedelta(days = 1)
	
	for delta in plotdays:
		filename = 'osisaf-average-last-' + str(delta) + '-days.png'
		fastSum(yesterday - timedelta(days = delta - 2), yesterday, filename)
		if putOnDropbox:
			uploadToDropbox([filename])
		os.remove(filename)
		
		filenameAntarctic = 'osisaf-antarctic-average-last-' + str(delta) + '-days.png'
		sumAntarctic(yesterday - timedelta(days = delta - 2), yesterday, filenameAntarctic)
		if putOnDropbox:
			uploadToDropbox([filenameAntarctic])
		os.remove(filenameAntarctic)

	if framExport:
		updateFramGraphs('fjl')
		updateFramGraphs('ateam')