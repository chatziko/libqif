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

def graphprivacy(dataset_name, gen_data): # adverror
	
	
	data_to_plot = []
	xticklabels = []
	for i in gen_data:
		content = []
		with open("generated_data_foursquare/"+i) as f:
			content = f.readlines()			
			print len(content)
		adverror = []
		for z in content:
			adverror.append(float(z))
		y = np.asarray(adverror)
		data_to_plot.append(y)

		l = i.split('-')
		if len(l) == 5:
			xticklabels.append('PL (Without labels)')
		else:
			xticklabels.append(l[1])



	fig = pyplot.figure(1,figsize=(21, 8))
	fig.suptitle(dataset_name+' Dataset', fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)

	bp = ax.boxplot(data_to_plot,widths=0.3,showmeans=True)
	ax.set_xticklabels(xticklabels)
	ax.set_xlabel('Semantic Labels')
	ax.set_ylabel('AdvError')
 	
 	pyplot.axhline(y=bp['boxes'][0].get_ydata()[2],linewidth=0.5)
	fig.savefig(dataset_name+'-foursquare.png', bbox_inches='tight')
	pyplot.show()



#################### Config ###########################


gen_data_tky = ['laplace-1.732868-tky.csv-foursquare-euclidean',
'laplace-Government Building-1.732868-tky.csv-foursquare-euclidean', 
'laplace-Airport-1.732868-tky.csv-foursquare-euclidean',
'laplace-Movie Theater-1.732868-tky.csv-foursquare-euclidean',
'laplace-Restaurant-1.732868-tky.csv-foursquare-euclidean',
'laplace-Indian Restaurant-1.732868-tky.csv-foursquare-euclidean',
'laplace-Dumpling Restaurant-1.732868-tky.csv-foursquare-euclidean',
'laplace-French Restaurant-1.732868-tky.csv-foursquare-euclidean',
'laplace-Dessert Shop-1.732868-tky.csv-foursquare-euclidean']

gen_data_nyc = ['laplace-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Airport-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Stadium-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Train Station-1.386294-nyc.csv-foursquare-euclidean',
'laplace-American Restaurant-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Italian Restaurant-1.386294-nyc.csv-foursquare-euclidean',
'laplace-French Restaurant-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Chinese Restaurant-1.386294-nyc.csv-foursquare-euclidean',
'laplace-Government Building-1.386294-nyc.csv-foursquare-euclidean']


#graphprivacy('Tokyo',gen_data_tky)
graphprivacy('New York',gen_data_nyc)

#graphlaplace()