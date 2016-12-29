#set size .5,1
set xtics (.1, .2, .3, .4, .5, .6, .7, .8, .9, 1)

set grid xtics
#set grid ytics
#set key left top

set lmargin 5
set rmargin 1
set bmargin 2

#set xrange [0:1]

set terminal pdf font 'Arial,22'
#set terminal epslatex size 9cm,7cm color colortext standalone header \
#   "\\newcommand{\\ft}[0]{\\footnotesize}"

#set parametric
#set trange [0:1.847]



#set xlabel "Probability p of first user"

#set ylabel "min-entropy leakage"
#set ytics .2
set yrange [1:100]
set output "crowds-repeated.pdf"
plot \
	"repeated.txt"		using 1:2	smooth frequency title "" with lines

