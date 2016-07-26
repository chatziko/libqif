import numpy as np
import matplotlib.pyplot as pyplot
import csv
def graphlaplace():

	dataset = "temp_foursquare/laplace-0.693147.csv"
	x = []
	for a in range(0,10000):
		x.append(a)
	x1  = np.asarray(x)
	with open(dataset) as csvfile:
		reader = csv.reader(csvfile)
		c = 0
		for row in reader:		
			y1 = np.asarray(row)
			pyplot.scatter(x1, y1, c = 'red')
			break			
			c = c+1
	pyplot.show()

def graphprivacy(): # adverror
	dataset1 = "generated_data_foursquare/laplace-1.732868-tky.csv-foursquare-euclidean"
	dataset2 = "generated_data_foursquare/laplace-Government Building-1.732868-tky.csv-foursquare-euclidean"
	dataset3 = "generated_data_foursquare/laplace-Airport-1.732868-tky.csv-foursquare-euclidean"
	dataset4 = "generated_data_foursquare/laplace-Movie Theater-1.732868-tky.csv-foursquare-euclidean"
	dataset6 = "generated_data_foursquare/laplace-Restaurant-1.732868-tky.csv-foursquare-euclidean"
	dataset7 = "generated_data_foursquare/laplace-Indian Restaurant-1.732868-tky.csv-foursquare-euclidean"
	dataset8 = "generated_data_foursquare/laplace-Dumpling Restaurant-1.732868-tky.csv-foursquare-euclidean"
	dataset9 = "generated_data_foursquare/laplace-French Restaurant-1.732868-tky.csv-foursquare-euclidean"
	dataset10 = "generated_data_foursquare/laplace-Dessert Shop-1.732868-tky.csv-foursquare-euclidean"

	dataset = [dataset1,dataset2,dataset3,dataset4,dataset6,dataset7,dataset8,dataset9,dataset10]

	data_to_plot = []

	for i in dataset:
		content = []
		with open(i) as f:
			content = f.readlines()
			print len(content)
		adverror = []
		for z in content:
			adverror.append(float(z))
		y = np.asarray(adverror)
		data_to_plot.append(y)


	fig = pyplot.figure(1,figsize=(14, 8))
	ax = fig.add_subplot(111)

	bp = ax.boxplot(data_to_plot,widths=0.4,showmeans=True)
	ax.set_xticklabels(['PL', 'Govt. Building','Airport','Movie Theater','Restaurant','Indian R','Dumpling R','French R','Dessert Shop'])
	
 	
 	pyplot.axhline(y=bp['boxes'][0].get_ydata()[2],linewidth=0.5)
	fig.savefig('tky_foursquare.png', bbox_inches='tight')
	pyplot.show()


graphprivacy()
#graphlaplace()