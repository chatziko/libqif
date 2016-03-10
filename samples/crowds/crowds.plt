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
set ytics .2
set yrange [-0.1:1]
set output "crowds-min-p.pdf"
plot \
	"data-min-p.txt"		using 1:2	smooth frequency title "phi = 0.0", \
	"data-min-p.txt"		using 1:3	smooth frequency title "phi = 0.5", \
	"data-min-p.txt"		using 1:4	smooth frequency title "phi = 1.0"

#set ylabel "tiger-leakage"
set yrange [-0.02:0.25]
set ytics .05
set output "crowds-tiger-p.pdf"
plot \
	"data-tiger-p.txt"		using 1:2	smooth frequency title "phi = 0.0", \
	"data-tiger-p.txt"		using 1:3	smooth frequency title "phi = 0.5", \
	"data-tiger-p.txt"		using 1:4	smooth frequency title "phi = 1.0"


#set xlabel "Probability phi of forwarding"

#set ylabel "min-entropy leakage"
set ytics .2
set yrange [-0.1:1]
set output "crowds-min-pf.pdf"
plot \
	"data-min-pf.txt"		using 1:2	smooth frequency title "p = 0.4", \
	"data-min-pf.txt"		using 1:3	smooth frequency title "p = 0.3", \
	"data-min-pf.txt"		using 1:4	smooth frequency title "p = 0.2"

#set ylabel "tiger-leakage"
set yrange [-0.02:0.25]
set ytics .05
set xtics (.1, .2, .3, .4, .5, .6, .75, .9, 1)
set output "crowds-tiger-pf.pdf"
plot \
	"data-tiger-pf.txt"		using 1:2	smooth frequency title "p = 0.4", \
	"data-tiger-pf.txt"		using 1:3	smooth frequency title "p = 0.3", \
	"data-tiger-pf.txt"		using 1:4	smooth frequency title "p = 0.2"

